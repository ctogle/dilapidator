
SceneJS.setConfigs({
	pluginPath: "scenejs_api/latest/plugins"
});

SceneJS.createScene({
	nodes: [
		{
			type: "cameras/pickFlyOrbit",
			yaw: -40,
			pitch: -20,
			zoom: 200,
			zoomSensitivity: 10.0,
			nodes: [
				{
					type: "rotate",
					x: 1,
					angle: -90,
					nodes: [
						{
    type: "texture",
    src: "world0/grass2.jpg",
    nodes: [
{
    type: "import/obj",
    src: "world0/model.5.mesh.obj",
},
{
    type: "import/obj",
    src: "world0/model.6.mesh.obj",
}
]
},
{
    type: "texture",
    src: "world0/orangeboxtex.png",
    nodes: [
{
    type: "import/obj",
    src: "world0/model.mesh.obj",
},
{
    type: "import/obj",
    src: "world0/model.2.mesh.obj",
},
{
    type: "import/obj",
    src: "world0/model.3.mesh.obj",
}
]
},
{
    type: "texture",
    src: "world0/concrete1.png",
    nodes: [
{
    type: "import/obj",
    src: "world0/model.4.mesh.obj",
}
]
}
					]
				}
			]
		}
	]
});