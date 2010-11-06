preprocess:
	@echo **********Creating mesh for each refinement.
	@$(foreach layers,4 8 16 32 64, \
		echo **********Installing directory for $(layers) layer simulation.;\
		install -d channel-$(layers); \
		cd channel-$(layers); \
		ln -s ../channel.flml .;\
		ln -s ../channel_tools.py .;\
		ln -s ../generate_mesh .;\
		echo **********Generating mesh for $(layers) layer simulation.;\
		./generate_mesh $(layers); \
		cd ..; \
	)

run:
	@echo **********Running simulation for each refinement.
	@$(foreach layers, 4 8 16 32 64, \
		cd channel-$(layers); \
		echo **********Running simulation for $(layers) layer mesh.;\
		echo ../../../bin/fluidity channel.flml; \
		../../../bin/fluidity channel.flml; \
		cd ..;\
	) 

postprocess:
	@echo **********Plotting convergence data.
	./plot_data 4 8 16 32 64

input: clean preprocess

clean:
	$(foreach layers,4 8 16 32 64, rm -r channel-$(layers);)
	rm -f *.stat
	rm -f *.d.* *.vtu
	rm -f *.msh *.geo
