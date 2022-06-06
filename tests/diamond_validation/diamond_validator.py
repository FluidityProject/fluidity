from itertools import chain
from xml.etree.ElementTree import ParseError

from diamond.schema import Schema


class DiamondValidator:
    def __init__(self):
        self.optionErrors = []

    def ValidateOptionsFiles(self, schemafile, testDir, extension, xmlRootNode):
        print(
            f"""\nValidating files against '{schemafile}' schema
Testing {extension} files with '{xmlRootNode}' root node"""
        )
        sch = Schema(schemafile)
        for filename in chain(*[list(testDir.rglob(f"*.{extension}"))]):
            try:
                optionsTree = sch.read(str(filename), testharness=True)
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
            if (
                any((lost_eles, added_eles, lost_attrs, added_attrs))
                or not optionsTree.valid
            ):
                self.optionErrors.append(filename.name)
                print(
                    f"""\n{"/".join(filename.parts[-2:])}: Errors
Lost elements: {lost_eles}
Added elements: {added_eles}
Lost attributes: {lost_attrs}
Added attributes: {added_attrs}
optionsTree.valid: {optionsTree.valid}"""
                )
