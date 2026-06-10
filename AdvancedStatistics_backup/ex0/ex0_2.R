# 1
# Initial attempt to load the data directly. 
# This fails to format correctly because the source is an HTML page, not a tabular file.
library(tibble)

url <- 'https://www.alltime-athletics.com/m_100ok.htm'
data <- read.csv(url, header=F)
as_tibble(data)

# 2
# Correct web approach. We extract the text inside the HTML <pre> tags 
# using 'rvest' and parse it as a fixed width file to perfectly capture the columns.
library(rvest)
library(tidyverse)
men100m_html <- read_html("http://www.alltime-athletics.com/m_100ok.htm")
men100m_html |> html_nodes(xpath = "//pre") |> html_text() -> men100m_list
men100m_tbl <- read_fwf(men100m_list)


# 3
# Data cleaning pipeline: we drop irrelevant bottom rows, convert text into 
# proper numeric/date formats, and rename variables for clarity.
men100m_tbl <- men100m_tbl |>
    slice(1:4848) |>
    mutate(X1 = as.integer(X1)) |>
    mutate(X2 = as.numeric(str_remove(X2, "A"))) |>
    mutate(X3 = as.numeric(str_remove(X3, "±"))) |>
    mutate(X9 = as.Date(X9, format = "%d.%m.%Y")) |>
    rename(
        Position = X1,
        Time = X2,
        Wind = X3,
        Athlete = X4,
        Nation = X5,
        Birth_Date = X6,
        Result_place = X7,
        Place = X8,
        Race_Date = X9
  )
print(men100m_tbl, n=10)
# NA values appear when there is no data

# 4
# a) Calculate the cumulative minimum of the times to track the historical 
# world record progression, then plot it as a line graph.
df_plot1 <- men100m_tbl |>
  arrange(Race_Date) |>
  mutate(AbsRecord = cummin(Time))
print(df_plot1[, c("Race_Date", "AbsRecord")], n=10)

plot(df_plot1$Race_Date, df_plot1$AbsRecord, 
     type = "l",                  # "l" stands for line
     col = "steelblue",         
     lwd = 1.5,                   # line width
     main = "Evolution of fastest time (Men)",
     xlab = "Years", 
     ylab = "Seconds (s)",
     panel.first = grid(col = "lightgray") # grid
)

# b) Identify the dominant nation per year by counting unique athletes, 
# and visualize the result using a dense, customized barplot.
top_nations_by_year <- men100m_tbl |> 
    drop_na(Race_Date) |>
    # Write each date with only the year
    mutate(Year = format(Race_Date, "%Y")) |>
    # .by groups by year and nation, and we create the new dataframe by adding
    # the number of athlete. We use n_dinstinct as an athlete can be present multiple times in a year 
    summarise(Unique_Runners = n_distinct(Athlete), .by = c(Year, Nation)) |>
    # for each year (column n=1) we select the maximum with respect to the number of unique runners
    slice_max(Unique_Runners, n = 1, by = Year) |>
    # we sort
    arrange(Year) |>
    # for plotting reason, we merge the year column with nation one in a new column
    mutate(Year_Country = paste(Year, Nation, sep = " - "))
top_nations_by_year

barplot(
    height = top_nations_by_year$Unique_Runners,
    names.arg = top_nations_by_year$Year_Country,
    col = "steelblue",
    main = "Number of Fastest Runners per year (Men)",
    xlab = "",
    ylab = "Counts",
    border = "black",
    space = 0.05,
    las = 2,                                          # rotate by 90 degrees
    cex.names = 0.5                                  # smaller the label text
)


# 5
# Same analysis with women dataset
url2 <- 'https://www.alltime-athletics.com/w_100ok.htm'

women100m_html <- read_html(url2)
women100m_html |> html_nodes(xpath = "//pre") |> html_text() -> women100m_list
women100m_tbl <- read_fwf(women100m_list)

women100m_tbl <- women100m_tbl |>
    slice(1:4848) |>
    mutate(X1 = as.integer(X1)) |>
    mutate(X2 = as.numeric(str_remove(X2, "A"))) |>
    mutate(X3 = as.numeric(str_remove(X3, "±"))) |>
    mutate(X9 = as.Date(X9, format = "%d.%m.%Y")) |>
    rename(
        Position = X1,
        Time = X2,
        Wind = X3,
        Athlete = X4,
        Nation = X5,
        Birth_Date = X6,
        Result_place = X7,
        Place = X8,
        Race_Date = X9
  )
print(women100m_tbl, n=10)

df_plot1 <- women100m_tbl |>
  arrange(Race_Date) |>
  mutate(AbsRecord = cummin(Time))
print(df_plot1[, c("Race_Date", "AbsRecord")], n=10)

plot(df_plot1$Race_Date, df_plot1$AbsRecord, 
     type = "l",                 
     col = "#f73b3b",         
     lwd = 2,                  
     main = "Evolution of fastest time (Women)",
     xlab = "Years", 
     ylab = "Seconds (s)",
     panel.first = grid(col = "lightgray")
)

top_nations_by_year <- women100m_tbl |> 
    drop_na(Race_Date) |>
    mutate(Year = format(Race_Date, "%Y")) |>
    summarise(Unique_Runners = n_distinct(Athlete), .by = c(Year, Nation)) |>
    slice_max(Unique_Runners, n = 1, by = Year) |>
    arrange(Year) |>
    mutate(Year_Country = paste(Year, Nation, sep = " - "))
top_nations_by_year

barplot(
    height = top_nations_by_year$Unique_Runners,
    names.arg = top_nations_by_year$Year_Country,
    col = "#f73b3b",
    main = "Number of Fastest Runners per year (Women)",
    xlab = "",
    ylab = "Counts",
    border = "black",
    space = 0.05,
    las = 2,
    cex.names = 0.5
)
