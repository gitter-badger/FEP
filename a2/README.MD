#Purpose
Reorders a PUMI based mesh using Reverse Cuthill McGee

You can obtain all of the test meshes from this other repository:

    $ git clone https://github.com/SCOREC/fep.git MESH_ROOT_DIR

The meshes we want will be in the folder **a2** of that project's master branch.

#Building
Just use make to build

    $ cd YOUR_PATH/FEP/a2
    $ make
    
There are three test meshes from availalbe in `MESH_ROOT_DIR/a2` that you cloned onto your local machine.
Run on a test mesh

    $ ./bin/a2 MESH_ROOT_DIR/a2/reorder_?.dmg MESH_ROOT_DIR/a2/reorder_?.smb
    
So to run on mesh "a" we will use:

    $ ./bin/a2 MESH_ROOT_DIR/a2/reorder_a.dmg MESH_ROOT_DIR/a2/reorder_a.smb
    
Note that we write the **.smb** filename **reorder_a.smb** even though the actual file name in `MESH_ROOT_DIR/a2` is something like **reorder_a0.smb**. This is due to how the PUMI library loads meshes in parallel. The number indicates which process should load which file. We are running only one process so there is only one part file: part 0.
