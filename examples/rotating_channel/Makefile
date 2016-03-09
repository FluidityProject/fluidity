include envcheck.mk

preprocess: envcheck
	@$(foreach layers,4 8 16 32 64, \
		echo **********Installing directory for $(layers) layer simulation.;\
		install -d channel-$(layers); \
		cd channel-$(layers); \
		ln -sf ../channel.flml .;\
		ln -sf ../channel_tools.py .;\
		ln -sf ../src/chan$(layers).msh ./channel.msh;\
		cd ..; \
	)

run: envcheck
	@echo **********Running simulation for each refinement.
	@$(foreach layers, 4 8 16 32 64, \
		cd channel-$(layers); \
		echo **********Running simulation for $(layers) layer mesh.;\
		echo ${BINPREFIX}fluidity -v 2 -l channel.flml; \
		${BINPREFIX}fluidity -v 2 -l channel.flml; \
		cd ..;\
	) 

postprocess: envcheck
	@echo **********Plotting convergence data.
	PYTHONPATH=${PWD}/../../python:$$PYTHONPATH ./plot_data 4 8 16 32 64

input: clean preprocess

clean:
	$(foreach layers,4 8 16 32 64, rm -rf channel-$(layers);)
	rm -f *.stat
	rm -f *.d.* *.vtu
