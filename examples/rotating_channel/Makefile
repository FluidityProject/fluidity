preprocess:
	@$(foreach layers,4 8 16 32 64, \
		echo **********Installing directory for $(layers) layer simulation.;\
		install -d channel-$(layers); \
		cd channel-$(layers); \
		ln -sf ../channel.flml .;\
		ln -sf ../channel_tools.py .;\
		ln -sf ../src/chan$(layers).msh ./channel.msh;\
		../../../bin/gmsh2triangle -2 ./channel.msh;\
		cd ..; \
	)

run:
	@echo **********Running simulation for each refinement.
	@$(foreach layers, 4 8 16 32 64, \
		cd channel-$(layers); \
		echo **********Running simulation for $(layers) layer mesh.;\
		echo ../../../bin/fluidity -v 2 -l channel.flml; \
		../../../bin/fluidity -v 2 -l channel.flml; \
		cd ..;\
	) 

postprocess:
	@echo **********Plotting convergence data.
	./plot_data 4 8 16 32 64

input: clean preprocess

clean:
	$(foreach layers,4 8 16 32 64, rm -rf channel-$(layers);)
	rm -f *.stat
	rm -f *.d.* *.vtu
