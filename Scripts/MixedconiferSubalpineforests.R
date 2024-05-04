
## ---------------------------
##
## Script name: 
##
## Author: Dr. Joan Dudney
##
## Date Created: 2024-04-23
##
## Copyright (c) Joan Dudney, 2024
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

theme_set(
  theme_bw(base_size = 15)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
)

summarize=dplyr::summarize
group_by=dplyr::group_by
select=dplyr::select

colors = c("#11302a","#036264", "#8f5774", "#e0829d","#dac4d0", "#6768ab", "#9389bd", "#cfcad1", "#2f2d3b", "#5e4d9b")


sierra <- data.frame(num = c(719, 127), loc = c("Mixed conifer", "Subalpine"))

sierra %>% 
  ggplot(aes(x = loc, y = num, fill = "blue", color = "darkgrey" ))+
  geom_bar(stat="identity", width = .5)+
  scale_fill_manual(values ="#036264")+
  scale_color_manual(values ="#cfcad1")+
  xlab("\nForest type in the Sierra Nevada")+
  ylab("# of publications")+
  guides(fill=F, col=F)+
  ylim (0, 1000)
  
