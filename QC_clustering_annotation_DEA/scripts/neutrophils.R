Idents(GSE141259) <- GSE141259$day



day7 <- subset(GSE141259,
                idents=7)

DimPlot(day7)

table(day10$cell.type, day10$group)

Idents(day7) <- day7$cell.type

DimPlot(day7)

neutro_day