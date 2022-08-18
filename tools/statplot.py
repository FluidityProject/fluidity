#!/usr/bin/env python3
import argparse
import warnings

import gi

gi.require_version("Gdk", "4.0")
gi.require_version("Gtk", "4.0")
import matplotlib.ticker as tck  # noqa: E402
import numpy as np  # noqa: E402
from gi.repository import Gdk, Gtk  # noqa: E402
from lxml import etree  # noqa: E402
from matplotlib.backends.backend_gtk4 import (  # noqa: E402
    NavigationToolbar2GTK4 as NavigationToolbar,
)
from matplotlib.backends.backend_gtk4agg import (  # noqa: E402
    FigureCanvasGTK4Agg as FigureCanvas,
)
from matplotlib.colors import TABLEAU_COLORS  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402


class Statplot(Gtk.ApplicationWindow):
    PLOT_COLOURS = list(TABLEAU_COLORS.keys())
    PLOT_MARKERS = Line2D.filled_markers[::-1]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_default_size(1600, 900)
        self.set_title("statplot")

        vertical_box = Gtk.Box(
            orientation=Gtk.Orientation.VERTICAL, homogeneous=False, spacing=8
        )
        self.set_child(vertical_box)

        self.entries_list, self.entry_values_list = self.extract_statfile_data()
        sorted_unique_entries = sorted(set().union(*self.entries_list))

        self.drop_down_x = self.create_drop_down(sorted_unique_entries, "ElapsedTime")
        self.drop_down_y = self.create_drop_down(sorted_unique_entries)

        horizontal_box = Gtk.Box(
            orientation=Gtk.Orientation.HORIZONTAL, homogeneous=True, spacing=8
        )
        horizontal_box.append(self.drop_down_x)
        horizontal_box.append(self.drop_down_y)

        self.fig = Figure(figsize=(14, 7), constrained_layout=True)
        self.ax = self.fig.add_subplot(1, 1, 1)

        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)

        vertical_box.append(horizontal_box)
        vertical_box.append(self.canvas)
        vertical_box.append(self.toolbar)

        self.add_event_controller_key("key-pressed", self.key_pressed, self.canvas)
        self.add_event_controller_key("key-pressed", self.key_pressed, self.toolbar)

        self.plot_type = "lines"
        self.update_plot_data()
        self.update_plot_axis("x", "linear")
        self.update_plot_axis("y", "linear")

        self.canvas.grab_focus()

    @staticmethod
    def add_event_controller_key(event, function, gtk_object):
        event_controller_key = Gtk.EventControllerKey.new()
        event_controller_key.connect(event, function)
        gtk_object.add_controller(event_controller_key)

    def create_drop_down(self, string_list, selected_item_string=None):
        string_list_items = "\n".ljust(11).join(
            [f"<item>{string}</item>" for string in string_list]
        )

        drop_down_ui_string = f"""<interface>
  <object class="GtkDropDown" id="drop-down">
    <property name="margin_start">100</property>
    <property name="margin_end">100</property>
    <property name="model">
      <object class="GtkStringList" id="string-list">
        <items>
          {string_list_items}
        </items>
      </object>
    </property>
    <property name="enable-search">true</property>
    <property name="expression">
      <lookup type="GtkStringObject" name="string"></lookup>
    </property>
  </object>
</interface>"""

        builder = Gtk.Builder.new_from_string(drop_down_ui_string, -1)
        drop_down = builder.get_object("drop-down")

        if selected_item_string is not None:
            drop_down.set_selected(string_list.index(selected_item_string))

        drop_down.connect("notify", self.drop_down_changed)

        return drop_down

    def drop_down_changed(self, drop_down, property):
        if property.name != "selected-item":
            return

        self.update_plot_data()
        self.update_plot_axis("x" if drop_down is self.drop_down_x else "y", "linear")

        self.canvas.grab_focus()

    def extract_statfile_data(self):
        entries_list, entry_values_list = [], []
        for statfile in self.statfile:
            entries, number_lines_header = self.gather_statfile_entries(statfile)
            entries_list.append(entries)
            entry_values_list.append(
                np.genfromtxt(statfile, skip_header=number_lines_header + 1)
            )
        return entries_list, entry_values_list

    @staticmethod
    def gather_statfile_entries(statfile):
        print(f"Reading {statfile}")

        entries = []
        with open(statfile, "r") as fid:
            for line_number, line in enumerate(fid):
                if line.startswith("</header>"):
                    break

                if not line.startswith("<field"):
                    continue

                field = etree.fromstring(line)

                if field.attrib["statistic"] == "value":
                    entries.append(field.attrib["name"])
                elif "components" in field.attrib and "material_phase" in field.attrib:
                    for component in range(int(field.attrib["components"])):
                        entries.append(
                            f"{field.attrib['material_phase']}%{field.attrib['name']}"
                            f"%{field.attrib['statistic']}%{component}"
                        )
                elif "components" in field.attrib:
                    for component in range(int(field.get("components"))):
                        entries.append(
                            f"{field.attrib['name']}"
                            f"%{field.attrib['statistic']}%{component}"
                        )
                elif "material_phase" in field.attrib:
                    entries.append(
                        f"{field.attrib['material_phase']}%{field.attrib['name']}"
                        f"%{field.attrib['statistic']}"
                    )
                else:
                    entries.append(
                        f"{field.attrib['name']}%{field.attrib['statistic']}"
                    )

        return entries, line_number

    def key_pressed(self, event, key_value, key_code, state):
        match key_value:
            case Gdk.KEY_l:
                if (
                    len(self.entry_values_list) == 1
                    and self.entry_values_list[0].ndim == 1
                ):
                    warnings.warn(
                        "Insufficient data available to turn on line display",
                        stacklevel=2,
                    )
                    return

                self.plot_type = "lines" if self.plot_type == "markers" else "markers"
                for artist_index, artist in enumerate(self.statplot):
                    if not artist.marker_only:
                        self.update_plot_type(artist, artist_index)
            case Gdk.KEY_q:
                self.destroy()
            case Gdk.KEY_r:
                self.entries_list, self.entry_values_list = self.extract_statfile_data()
                sorted_unique_entries = sorted(set().union(*self.entries_list))

                current_x_entry = self.drop_down_x.get_selected_item().get_string()
                current_y_entry = self.drop_down_y.get_selected_item().get_string()

                self.drop_down_x.set_model(Gtk.StringList.new(sorted_unique_entries))
                self.drop_down_y.set_model(Gtk.StringList.new(sorted_unique_entries))

                self.drop_down_x.set_selected(
                    sorted_unique_entries.index(current_x_entry)
                )
                self.drop_down_y.set_selected(
                    sorted_unique_entries.index(current_y_entry)
                )

                self.update_plot_data()
            case Gdk.KEY_x | Gdk.KEY_y:
                if (
                    len(self.entry_values_list) == 1
                    and self.entry_values_list[0].ndim == 1
                ):
                    warnings.warn(
                        "Insufficient data available to turn on logarithmic scale",
                        stacklevel=2,
                    )
                    return

                if key_value == Gdk.KEY_x:
                    self.get_scale = self.ax.get_xscale
                    axis_selected = "x"
                    current_data = np.hstack(self.plot_data_x)
                else:
                    self.get_scale = self.ax.get_yscale
                    axis_selected = "y"
                    current_data = np.hstack(self.plot_data_y)

                match self.get_scale():
                    case "linear":
                        data_min = min(current_data)
                        data_max = max(current_data)

                        if data_min == data_max:
                            warnings.warn(
                                "Change to logarithmic scale denied: the selected "
                                f"variable for the {axis_selected} axis is a constant.",
                                stacklevel=2,
                            )
                            return
                        elif data_min < 0:
                            scale = "symlog"
                        else:
                            scale = "log"
                    case "log" | "symlog":
                        scale = "linear"
                    case _:
                        raise SystemExit(f"Unknown axis scale: {self.get_scale()}")

                self.update_plot_axis(axis_selected, scale)
            case Gdk.KEY_Escape:
                self.canvas.grab_focus()

    def update_plot_axis(self, axis, scale):
        match axis:
            case "x":
                current_axis = self.ax.xaxis
                current_data = np.hstack(self.plot_data_x)
                drop_down = self.drop_down_x
                set_scale = self.ax.set_xscale
                set_label = self.ax.set_xlabel
            case "y":
                current_axis = self.ax.yaxis
                current_data = np.hstack(self.plot_data_y)
                drop_down = self.drop_down_y
                set_scale = self.ax.set_yscale
                set_label = self.ax.set_ylabel
            case _:
                raise SystemExit(f"Unknown axis: {axis}")

        match scale:
            case "linear":
                set_scale(scale)
                self.ax.ticklabel_format(
                    style="sci", axis=axis, scilimits=(0, 0), useMathText=True
                )
                current_axis.set_minor_locator(tck.AutoMinorLocator())
            case "log":
                set_scale(scale)
            case "symlog":
                data_min = min(abs(current_data[current_data != 0]))
                data_max = max(abs(current_data))
                axis_range = np.log10(data_max / data_min)

                set_scale(
                    scale,
                    base=10,
                    linthresh=data_min * 10 ** (axis_range / 2),
                    subs=np.arange(2, 10),
                    linscale=0.5,
                )
            case _:
                raise SystemExit(f"Unknown axis scale: {scale}")

        set_label(
            drop_down.get_selected_item().get_string(), fontweight="bold", fontsize=18
        )

        self.ax.tick_params(axis=axis, which="major", length=8, width=3, labelsize=16)
        self.ax.tick_params(
            axis=axis, which="minor", length=4, width=2, labelsize=12, colors="crimson"
        )

        current_axis.get_offset_text().set(
            fontsize=14, fontweight="bold", color="black"
        )

        self.ax.relim()
        self.ax.autoscale()

        self.fig.canvas.draw_idle()

    def update_plot_data(self):
        try:
            self.plot_data_x.clear()
            self.plot_data_y.clear()
        except AttributeError:
            self.plot_data_x, self.plot_data_y, self.statplot = [], [], []

        artist_index = 0
        current_artist_number = len(self.statplot)

        for entries, entry_values in zip(self.entries_list, self.entry_values_list):
            try:
                field_entry_index_x = entries.index(
                    self.drop_down_x.get_selected_item().get_string()
                )
                field_entry_index_y = entries.index(
                    self.drop_down_y.get_selected_item().get_string()
                )

                self.plot_data_x.append(
                    np.atleast_1d(entry_values[..., field_entry_index_x])
                )
                self.plot_data_y.append(
                    np.atleast_1d(entry_values[..., field_entry_index_y])
                )

                single_point = self.plot_data_x[-1].size == 1

                if artist_index < current_artist_number:
                    self.statplot[artist_index].set_xdata(self.plot_data_x[-1])
                    self.statplot[artist_index].set_ydata(self.plot_data_y[-1])

                    if (
                        self.statplot[artist_index].marker_only ^ single_point
                    ) and self.plot_type == "lines":
                        self.update_plot_type(self.statplot[artist_index], artist_index)
                else:
                    (artist,) = self.ax.plot(
                        self.plot_data_x[-1],
                        self.plot_data_y[-1],
                        linestyle="solid"
                        if single_point or self.plot_type == "markers"
                        else "None",
                    )
                    self.update_plot_type(artist, artist_index)

                    self.statplot.append(artist)

                self.statplot[artist_index].marker_only = single_point

                artist_index += 1

            except ValueError:
                continue

        for artist in self.statplot[artist_index:]:
            artist.remove()
        del self.statplot[artist_index:]

        self.toolbar.update()
        self.fig.canvas.draw_idle()

    def update_plot_type(self, artist, artist_index):
        if artist.get_linestyle() == "None":
            artist.set_color(self.PLOT_COLOURS[artist_index % len(self.PLOT_COLOURS)])
            artist.set_linestyle("solid")
            artist.set_linewidth(3)
            artist.set_marker("None")
        else:
            artist.set_linestyle("None")
            artist.set_marker(self.PLOT_MARKERS[artist_index % len(self.PLOT_MARKERS)])
            artist.set_markeredgecolor("black")
            artist.set_markeredgewidth(0.5)
            artist.set_markerfacecolor(
                self.PLOT_COLOURS[artist_index % len(self.PLOT_COLOURS)]
            )
            artist.set_markersize(10)

        self.fig.canvas.draw_idle()


class StatApp(Gtk.Application):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.connect("activate", self.on_activate)

    def on_activate(self, app):
        self.win = Statplot(application=app)
        self.win.show()


parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description="""statplot, a GUI to visualise Fluidity .stat outputs

When the window is displayed, press the following keys for additional options:
- x, y -> Change the scale of the x, y axis (linear or logarithmic).
- l    -> Change how data is represented (solid line or markers).
- r    -> Reload the provided .stat file(s).
- q    -> Close the GUI.""",
)
parser.add_argument("statfile", nargs="+", help="Path(s) to the .stat file(s)")
args = parser.parse_args()

Statplot.statfile = args.statfile
app = StatApp()
app.run()
