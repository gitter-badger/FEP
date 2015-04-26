#Purpose
Reorders a PUMI based mesh using reverse cuthill McGee

You can obtain all of the test meshes from this repository:

    $ git clone https://github.com/SCOREC/fep.git MESH_ROOT_DIR

The meshes we want will be in the folder **a2** of the master branch.

#Building
To build and run on a test mesh

    $ make
    $ ./bin/a2 MESH_ROOT_DIR/a2/reorder_a.dmg MESH_ROOT_DIR/a2/reorder_a.smb
