dt_A = 1.0 
dt_B = 0.5
dt_C = 0.25
dt_D = 0.125
dt_E = 0.0625
source_swml = wave_periodic_source.swml

input: clean copy_tests; cd src; make input

run: copy_tests run_dt run_dt_dx

run_dt:
	shallow_water wave_dt_A.swml; shallow_water wave_dt_B.swml; shallow_water wave_dt_C.swml; shallow_water wave_dt_D.swml

run_dt_dx:
	shallow_water wave_dt_dx_A.swml; shallow_water wave_dt_dx_B.swml; shallow_water wave_dt_dx_C.swml; shallow_water wave_dt_dx_D.swml

run_parallel:
	(../../bin/shallow_water wave_dt_A.swml &); (../../bin/shallow_water wave_dt_B.swml &); (../../bin/shallow_water wave_dt_C.swml &); (../../bin/shallow_water wave_dt_D.swml &); (../../bin/shallow_water wave_dt_dx_A.swml &); (../../bin/shallow_water wave_dt_dx_B.swml &); (../../bin/shallow_water wave_dt_dx_C.swml &); (../../bin/shallow_water wave_dt_dx_D.swml &); wait

copy_tests:
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_A)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/wave/wave_dt_A/' | sed 's/src\/mesh/src\/mesh_D/' > wave_dt_A.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_B)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/wave/wave_dt_B/' | sed 's/src\/mesh/src\/mesh_D/' > wave_dt_B.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_C)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/wave/wave_dt_C/' | sed 's/src\/mesh/src\/mesh_D/' > wave_dt_C.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_D)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/wave/wave_dt_D/' | sed 's/src\/mesh/src\/mesh_D/' > wave_dt_D.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_A)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/wave/wave_dt_dx_A/' | sed 's/src\/mesh/src\/mesh_A/' > wave_dt_dx_A.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_B)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/wave/wave_dt_dx_B/' | sed 's/src\/mesh/src\/mesh_B/' > wave_dt_dx_B.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_C)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/wave/wave_dt_dx_C/' | sed 's/src\/mesh/src\/mesh_C/' > wave_dt_dx_C.swml
	sed -n '1h;1!H;$${;g;s/<timestep>.*<\/timestep>/<timestep>\n<real_value rank="0">$(dt_D)<\/real_value>\n<\/timestep>/g;p;}' $(source_swml) | sed 's/wave/wave_dt_dx_D/' | sed 's/src\/mesh/src\/mesh_D/' > wave_dt_dx_D.swml
clean:
	rm -f *.vtu *.stat wave_dt_*.swml src/mesh*.msh 
