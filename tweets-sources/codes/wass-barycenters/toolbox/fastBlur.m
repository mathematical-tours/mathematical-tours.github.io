function result = fastBlur(image, filterWidth)

result = fastBlur1d(fastBlur1d(image,filterWidth)', filterWidth)';