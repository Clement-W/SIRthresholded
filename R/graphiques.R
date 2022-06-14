
# p le nombre de variables dans X
barplot_model_size = function(listevarHard,listevarSoft) {

    modelsizeHard = rep(0, length(listevarHard))
    for (i in 1:length(listevarHard)) {
        modelsizeHard[i] = length(listevarHard[[i]])
    }
    
    modelsizeSoft = rep(0, length(listevarSoft))
    for (i in 1:length(listevarSoft)) {
        modelsizeSoft[i] = length(listevarSoft[[i]])
    }
    
    #variables manquantes dans hard:
    varManquantesHard = modelsizeSoft[which(!modelsizeSoft %in% modelsizeHard)]
    modelsizeHard = c(modelsizeHard,varManquantesHard)
    
    #variables manquantes dans soft:
    varManquantesSoft = modelsizeHard[which(!modelsizeHard %in% modelsizeSoft)]
    modelsizeSoft = c(modelsizeSoft,varManquantesSoft)

    tableHard = table(modelsizeHard) / length(listevarHard)
    tableSoft = table(modelsizeSoft) / length(listevarSoft)
    
    tableHard[which(names(tableHard) %in% paste(varManquantesHard))] = NA
    tableSoft[which(names(tableSoft) %in% paste(varManquantesSoft))] = NA
    
    tables = rbind(tableHard,tableSoft)*100   

    par(mar=c(5, 4, 6, 2), xpd=TRUE)
    barplot(tables, xlab = "Size of the reduced model", ylab = "percent",beside = T,col=c("steelblue3","tan2"))
    legend("topright",legend = c("HT-SIR","ST-SIR"),pch=16,col=c("steelblue3","tan2"),cex=1.1,inset = c(0,-0.2))

}

barplot_occurrences = function(listevarHard,listevarSoft) {
    nbSim = length(listevarHard)
    
    listevarHard = unlist(listevarHard)
    listevarSoft = unlist(listevarSoft)
    
    
    #variables manquantes dans hard:
    varManquantesHard = listevarSoft[which(!listevarSoft %in% listevarHard)]
    listevarHard = c(listevarHard,varManquantesHard)
    
    #variables manquantes dans soft:
    varManquantesSoft = listevarSoft[which(!listevarHard %in% listevarSoft)]
    listevarSoft = c(listevarSoft,varManquantesSoft)

    tableHard = table(listevarHard) / nbSim
    valsHard <- as.numeric(gsub("X","", names(tableHard)))
    tableHard = tableHard[order(valsHard)]
    
    tableSoft = table(listevarSoft) / nbSim
    valsSoft <- as.numeric(gsub("X","", names(tableSoft)))
    tableSoft = tableSoft[order(valsSoft)]
    
    tableHard[which(names(tableHard) %in% paste(varManquantesHard))] = NA
    tableSoft[which(names(tableSoft) %in% paste(varManquantesSoft))] = NA

    tables = rbind(tableHard,tableSoft)*100 

    par(mar=c(5, 4, 6, 2), xpd=TRUE)
    barplot(tables, xlab = "selected variables", ylab = "percent",beside = T,col=c("steelblue3","tan2"))
    legend("topright",legend = c("HT-SIR","ST-SIR"),pch=16,col=c("steelblue3","tan2"),cex=1.1,inset = c(0,-0.2))
}

boxplot_cos2 = function(betacos2){
    boxplot(list(betacos2[, 1], betacos2[, 2], betacos2[, 3], betacos2[, 4]), names = c("SIR", "HT-SIR", "ST-SIR", "ST-SIR-p*"),ylab="cos²",col=c("gray","steelblue3","tan2","darkolivegreen3"))
    #title(paste("cos^2 entre beta théorique et beta estimé selon les méthodes\n", "n=", nn, ", p=", pp, ", p1=", pp1, ", ratio=", RRATIO, "link fct=", fct, "CorrelationX=", correlationX))
}

nbVarSelectionnees = function(listevarHard,listevarSoft){
    nbVarHard = rep(0,length(listevarHard))
    nbVarSoft = rep(0,length(listevarSoft))
    for (i in 1:length(listevarHard)){
        nbVarSoft[i] = length(listevarSoft[[i]])
        nbVarHard[i] = length(listevarHard[[i]])
    }
    boxplot(nbVarHard, nbVarSoft, names = c("HT-SIR", "ST-SIR"))
}
