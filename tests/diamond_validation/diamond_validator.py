from diamond.schema import Schema
from itertools import chain
from xml.etree.ElementTree import ParseError


class DiamondValidator:
    def __init__(self):
        self.optionErrors = []

    def ValidateOptionsFiles(self, schemafile, testDir, extension,
                             xmlRootNode):
        print(f"""\nValidating files against '{schemafile}' schema
Testing {extension} files with '{xmlRootNode}' root node""")
        sch = Schema(schemafile)
        for filename in chain(*[list(testDir.rglob(f"*/*.{extension}"))]):
            try:
                optionsTree = sch.read(str(filename))
            except ParseError:
                self.optionErrors.append(filename)
                print(f"\n{filename.name}: Parsing error")
                continue
            else:
                if optionsTree is None:
                    self.optionErrors.append(filename.name)
                    print(f"\n{filename.name}: Syntax error")
                    continue
            lost_eles, added_eles, lost_attrs, added_attrs = sch.read_errors()
            if (len(lost_eles + added_eles + lost_attrs + added_attrs) > 0
                    or not optionsTree.valid):
                self.optionErrors.append(filename.name)
                print(f"""\n{filename.name}: Errors
Lost elements: {lost_eles}
Added elements: {added_eles}
Lost attributes: {lost_attrs}
Added attributes: {added_attrs}
optionsTree.valid: {optionsTree.valid}""")
