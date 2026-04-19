# 1
# We define the vectors corresponding to the table's columns and then the dataframe.
names = c("Loch Ness", "Loch Lomond", "Loch Morar", "Loch Tay", "Loch Awe", "Loch Maree", "Loch Ericht", "Loch Lochy", "Loch Rannoch", "Loch Shiel", "Loch Katrine", "Loch Arkaig", "Loch Shin")
volumes = c(7.45, 2.6, 2.3, 1.6, 1.2, 1.09, 1.08, 1.07, 0.97, 0.79, 0.77, 0.75, 0.35)
areas = c(56, 71, 27, 26.4, 39, 28.6, 18.6, 16, 19, 19.5, 12.4, 16, 22.5)
lengths = c(39, 36, 18.8, 23, 41, 20, 23, 16, 15.7, 28, 12.9, 19.3, 27.8)
max_depths = c(230, 190, 310, 150, 94, 114, 156, 162, 134, 128, 151, 109, 49)
mean_depths = c(132, 37, 87, 60.6, 32, 38, 57.6, 70, 51, 40, 43.4, 46.5, 15.5)

scottish.lakes <- data.frame(names, volumes, areas, lengths, max_depths, mean_depths)
scottish.lakes

# 2
# We determine the lake with largest/smallest volume/surface using
# which.max and which.min that select the index of maximum/minimum value of a vector
# Then we get the lake's name by accessing the df attribute.
largest_volume <- scottish.lakes$names[which.max(scottish.lakes$volumes)]
largest_volume
smallest_volume <- scottish.lakes$names[which.min(scottish.lakes$volumes)]
smallest_volume

largest_surface <- scottish.lakes$names[which.max(scottish.lakes$areas)]
largest_surface
smallest_surface <- scottish.lakes$names[which.min(scottish.lakes$areas)]
smallest_surface

# 3
# We sort the dataframe by using order, which sorts the selected vector and returns
# the element indexes. Then, since the elements are in ascending order, we tail the last two
# elements to get the two largest lakes.
scottish.lakes_sorted <- scottish.lakes[order(scottish.lakes$areas), ]
scottish.lakes_sorted
tail(scottish.lakes_sorted, 2)$names

# 4
# The total area is simply the sum of the areas.
sum(scottish.lakes$areas)

