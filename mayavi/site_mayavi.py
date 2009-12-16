from enthought.mayavi.core.registry import registry
from enthought.mayavi.core.pipeline_info import PipelineInfo
from enthought.mayavi.core.metadata import FilterMetadata
from enthought.mayavi.core.metadata import SourceMetadata

# Metadata for the new filters we want to add
boundary_marker_editor = FilterMetadata(
    id         = "BoundaryMarkerEditor",
    menu_name  = "BoundaryMarkerEditor",
    factory    = 'mayavi_amcg.filters.boundary_marker_editor.BoundaryMarkerEditor'
)

field_operations = FilterMetadata(
    id         = "FieldOperations",
    menu_name  = "FieldOperations",
    factory    = 'mayavi_amcg.filters.field_operations.FieldOperations'
)

projection_and_depth_stretch = FilterMetadata(
    id         = "ProjectionAndDepthStretch",
    menu_name  = "ProjectionAndDepthStretch",
    factory    = 'mayavi_amcg.filters.projection_and_depth_stretch.ProjectionAndDepthStretch'
)

mesh_diagnostics = FilterMetadata(
    id         = "MeshDiagnostics",
    menu_name  = "MeshDiagnostics",
    factory    = 'mayavi_amcg.filters.mesh_diagnostics.MeshDiagnostics'
)

tensor_eigenvectors_eigenvalues = FilterMetadata(
    id         = "TensorEigenvectorsEigenvalues",
    menu_name  = "TensorEigenvectorsEigenvalues",
    factory    = 'mayavi_amcg.filters.tensor_eigenvectors_eigenvalues.TensorEigenvectorsEigenvalues'
)

# Register the filters with the mayavi registry
registry.filters.append(boundary_marker_editor)
registry.filters.append(field_operations)
registry.filters.append(projection_and_depth_stretch)
registry.filters.append(mesh_diagnostics)
registry.filters.append(tensor_eigenvectors_eigenvalues)

# Metadata for the new source we want to add
triangle_reader_info = SourceMetadata(
    id          = "TriangleReader",
    class_name  = 'mayavi_amcg.triangle_reader.TriangleReader',
    tooltip     = "Load Triangle files",
    desc        = "Load Triangle files",
    help        = "Load Triangle files",
    menu_name   = "&Triangle files",
    extensions  = ['face','edge','ele'],
    wildcard    = 'Triangle files (*.face)|*.face|Triangle files (*.edge)|*.edge|Triangle files (*.ele)|*.ele',
    output_info = PipelineInfo(datasets=['unstructured_grid'],
                               attribute_types=['any'],
                               attributes=['any'])
)

# Register the source with the mayavi registry
registry.sources.append(triangle_reader_info)


if __name__ == '__main__':
    import sys
    print "*"*80
    print "ERROR: This script isn't supposed to be executed."
    print __doc__
    print "*"*80

    from enthought.util.home_directory import get_home_directory
    print "Your .mayavi2 directory should be in %s"%get_home_directory()
    print "*"*80
    sys.exit(1)
