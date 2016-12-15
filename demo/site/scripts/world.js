// Using the "import/obj" node to import a raptor mesh from .OBJ format
// Internally, the node uses the K3D library for parsing - http://k3d.ivank.net/

// Point SceneJS to the bundled plugins
SceneJS.setConfigs({
    pluginPath: "scenejs/api/latest/plugins"
});

// Create a scene
SceneJS.createScene({
        nodes: [

            // Mouse-orbited camera,implemented by plugin at
            // http://scenejs.org/api/latest/plugins/node/cameras/orbit.js
            {
                type: "cameras/orbit",
                yaw: -40,
                pitch: -20,
                zoom: 200,
                zoomSensitivity: 20.0,

                nodes: [

                    // Move the raptor a bit to centre it
                    {
                        type: "translate", y: -30, z: -20,
                        nodes: [

                            // Texture our raptor - the "import.obj" node merely imports meshes, so
                            // it's our job to apply any materials and textures. This tends to be OK
                            // because for efficiency we would usually have one big baked texture
                            // like this one.
                            //
                            // Texture is loaded from:
                            // http://scenejs.org/examples/models/obj/raptor.jpg

                            {
                                type: "texture",
                                src: "world0/orangeboxtex.png",

                                nodes: [
                                    // Import the .OBJ mesh
                                    //
                                    // This node is implemented by plugin at:
                                    // http://scenejs.org/api/latest/plugins/node/import/obj.js
                                    //
                                    // The OBJ file is loaded from:
                                    // http://scenejs.org/examples/models/obj/raptor.obj
                                    {
                                        type: "import/obj",
                                        src: "world0/model.mesh.obj",
                                    }
                                ]
                            }
                        ]
                    }
                ]
            }
        ]
    }
);
