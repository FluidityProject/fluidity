dt_A = 0.1
dt_B = 0.05
dt_C = 0.025
dt_D = 0.0125
dt_E = 0.006125
source_swml = steady_periodic.swml

input: clean copy_tests
	cd src; make

run: copy_tests run_dt run_dt_dx

run_dt:
	shallow_water steady_dt_A.swml; shallow_water steady_dt_B.swml; shallow_water steady_dt_C.swml; shallow_water steady_dt_D.swml; shallow_water steady_dt_E.swml

run_dt_dx:
	shallow_water steady_dt_dx_A.swml; shallow_water steady_dt_dx_B.swml; shallow_water steady_dt_dx_C.swml; shallow_water steady_dt_dx_D.swml; shallow_water steady_dt_dx_E.swml

copy_tests:
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_A)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/steady_periodic/steady_dt_A/' | sed 's/src\/mesh/src\/mesh_E/' > steady_dt_A.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_B)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/steady_periodic/steady_dt_B/' | sed 's/src\/mesh/src\/mesh_E/' > steady_dt_B.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_C)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/steady_periodic/steady_dt_C/' | sed 's/src\/mesh/src\/mesh_E/' > steady_dt_C.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_D)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/steady_periodic/steady_dt_D/' | sed 's/src\/mesh/src\/mesh_E/' > steady_dt_D.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_E)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/steady_periodic/steady_dt_E/' | sed 's/src\/mesh/src\/mesh_E/' > steady_dt_E.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_A)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/steady_periodic/steady_dt_dx_A/' | sed 's/src\/mesh/src\/mesh_A/' > steady_dt_dx_A.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_B)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/steady_periodic/steady_dt_dx_B/' | sed 's/src\/mesh/src\/mesh_B/' > steady_dt_dx_B.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_C)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/steady_periodic/steady_dt_dx_C/' | sed 's/src\/mesh/src\/mesh_C/' > steady_dt_dx_C.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_D)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/steady_periodic/steady_dt_dx_D/' | sed 's/src\/mesh/src\/mesh_D/' > steady_dt_dx_D.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_E)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/steady_periodic/steady_dt_dx_E/' | sed 's/src\/mesh/src\/mesh_E/' > steady_dt_dx_E.swml
clean:
	rm -f *.vtu *.stat steady_dt_*.swml
	cd src; rm mesh*
