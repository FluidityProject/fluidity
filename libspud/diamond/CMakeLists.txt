configure_file(setup.py.in setup.py)
execute_process(COMMAND ${Python3_EXECUTABLE} setup.py build)
install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} setup.py install)")

# https://stackoverflow.com/questions/29232614/build-a-python-package-with-setup-py-in-cmake
