# Spud

Spud is a generic system for defining, writing and processing options files for scientific computer models.

The interfaces to scientific computer models are frequently primitive, under-documented and ad-hoc text files. This makes using and developing the model in question difficult and error-prone.

With Spud, the model developer need only write a rules file (schema) which defines the options which the model takes and the relationship between them. The Spud component Diamond then provides an automatically generated graphical user interface which guides the user and validates the user's input against the schema. Diamond writes out an xml options file for use in Spud.

The developer then uses libspud to read the options file into the model. Libspud can read any valid options file without further code modifications and makes the options available at any point in the model code at which they are required.

Spud further provides the facility for the schema to be self-documenting and Diamond presents this documentation to the model user in a context-sensitive manner.

## Documentation

For installation and usage instructions, please see doc/spud_manual.pdf.

For an example of spud and diamond usage, see the examples directory.
