# script.R
# Ashley Mae Conard
# Last Mod. 6 June 2020

app <- shinyApp(ui = ui, server = server)
runApp(app, host = "0.0.0.0", port = 80, launch.browser = FALSE)
