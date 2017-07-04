
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
    src: "../world/grass2.jpg",
    nodes: [
{
    type: "import/obj",
    src: "../world/model.mesh.obj",
}
]
}
					]
				}
			]
		}
	]
});