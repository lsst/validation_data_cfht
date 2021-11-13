# Skymap config for this CFHT data; one tract centered on the observations.
config.skyMap.name = "discrete"

config.skyMap['discrete'].projection='TAN'

# dimensions of inner region of patches (x,y pixels)
config.skyMap['discrete'].patchInnerDimensions=[4000, 4000]

# nominal pixel scale (arcsec/pixel)
config.skyMap['discrete'].pixelScale=0.263

# center of this data
config.skyMap['discrete'].raList=[214.856821]
config.skyMap['discrete'].decList=[52.662694]
config.skyMap['discrete'].radiusList=[0.6]
