DATA <- read.table("/home/itamblyn/fortran/error/data.dat", header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
linear_data <- lm(formula = V ~ I, data = DATA)          
summary(linear_data)

ERROR <- read.table("/home/itamblyn/fortran/error/error.dat", header=TRUE, sep="", na.strings="NA", dec=".", strip.white=TRUE)
linear_error <- lm(formula = V ~ I, data = ERROR)          
summary(linear_error)
