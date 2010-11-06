preprocess:
	@echo **********Creating mesh for each refinement N = 1/16 -- 1/128.
	@$(foreach NN,16 32 64 128, \
		echo **********Installing directory for N = $(NN) simulation.;\
		install -d driven_cavity-$(NN); \
		cd driven_cavity-$(NN); \
		ln -s ../driven_cavity.flml .;\
		ln -s ../driven_cavity.py .;\
		ln -s ../generate_mesh .;\
		echo **********Generating mesh for N = $(NN) simulation.;\
		./generate_mesh $(NN); \
		cd ..; \
	)

run:
	@echo **********Running simulation for each refinement.
	@$(foreach NN,16 32 64 128, \
		cd driven_cavity-$(NN); \
		echo **********Running simulation for N = $(NN) mesh.;\
		echo ../../../bin/fluidity -v2 -l driven_cavity.flml; \
		../../../bin/fluidity -v2 -l driven_cavity.flml; \
		cd ..;\
	) 

postprocess:
	@echo **********Plotting convergence data.
	./plot_data 16 32 64 128

input: clean preprocess

clean:
	$(foreach NN,16 32 64 128, rm -rf driven_cavity-$(NN);)
	rm -f *.stat
	rm -f *.d.* *.vtu
	rm -f *.msh *.geo
