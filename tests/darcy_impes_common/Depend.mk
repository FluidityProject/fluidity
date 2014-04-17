
define import_python_code
@echo Importing $1
@echo '#!/usr/bin/env python' > temp
@echo '' >> temp
@echo '# ************ THIS FILE IS A COPY THAT WILL BE DISCARDED. ************' >> temp
@echo '#   ' >> temp
@echo '#   The copy has been made as an alternative to using Python''s import '
@echo '#   function, as the latter doesn''t work well in the Fluidity '
@echo '#   test harness.  If you want to edit the file, consider editing '
@echo '#   ../darcy_impes_common/$1 instead.' >> temp
@echo '#   ' >> temp
@echo '# *********************************************************************' >> temp
@echo '' >> temp
@echo '' >> temp
@cat ../darcy_impes_common/$1 >> temp
@mv temp $1;
endef
