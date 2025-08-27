music_basis = function(x, non.zero = TRUE, markers = NULL, clusters, samples, select.ct = NULL, cell_size = NULL, ct.cov = FALSE, verbose = TRUE){
  if(!is.null(select.ct)){
    x = x[, x@colData[, clusters] %in% select.ct]
  }
  if(non.zero){  ## eliminate non expressed genes
    x <- x[rowSums(counts(x))>0, ]
  }
  
  clusters <- as.character(colData(x)[, clusters])
  samples <- as.character(colData(x)[, samples])
  
  M.theta <- sapply(unique(clusters), function(ct){
    my.rowMeans(sapply(unique(samples), function(sid){
      y = counts(x)[,clusters %in% ct & samples %in% sid]
      if(is.null(dim(y))){
        y <- as.vector(y)
        return(y/sum(y))
      }else{
        return(rowSums(y)/sum(y))
      }
    }), na.rm = TRUE)
  })
  if(verbose){message("Creating Relative Abudance Matrix...")}
  if(ct.cov){
    nGenes = nrow(x);
    n.ct = length(unique(clusters));
    nSubs = length(unique(samples))
    
    Theta <- sapply(unique(clusters), function(ct){
      sapply(unique(samples), function(sid){
        y = counts(x)[,clusters %in% ct & samples %in% sid]
        if(is.null(dim(y))){
          return(y/sum(y))
        }else{
          return( rowSums(y)/sum(y) )
        }
      })
    })
    if(!is.null(select.ct)){
      m.ct = match(select.ct, colnames(Theta))
      Theta = Theta[, m.ct]
    }
    
    Sigma.ct = sapply(1:nGenes, function(g){
      sigma.temp = Theta[nGenes*(0:(nSubs - 1)) + g, ];
      Cov.temp = cov(sigma.temp)
      Cov.temp1 = cov(sigma.temp[rowSums(is.na(Theta[nGenes*(0:(nSubs - 1)) + 1, ])) == 0, ])
      Cov.temp[which(colSums(is.na(sigma.temp))>0), ] = Cov.temp1[which(colSums(is.na(sigma.temp))>0), ]
      Cov.temp[, which(colSums(is.na(sigma.temp))>0)] = Cov.temp1[, which(colSums(is.na(sigma.temp))>0)]
      return(Cov.temp)
    })
    colnames(Sigma.ct) = rownames(x);
    
    if (!is.null(markers)){
      ids <- intersect(unlist(markers), rownames(x))
      m.ids = match(ids, rownames(x))
      Sigma.ct <- Sigma.ct[ , m.ids]
    }
    if(verbose){message("Creating Covariance Matrix...")}
  }else{
    Sigma <- sapply(unique(clusters), function(ct){
      apply(sapply(unique(samples), function(sid){
        y = counts(x)[,clusters %in% ct & samples %in% sid]
        if(is.null(dim(y))){
          y <- as.vector(y)
          vec <- y/sum(y)
          names(vec) <- rownames(x)
          return(vec)
        }else{
          return(rowSums(y)/sum(y))
        }
      }), 1, var, na.rm = TRUE)
    })
    if(!is.null(select.ct)){
      m.ct = match(select.ct, colnames(Sigma))
      Sigma = Sigma[, m.ct]
    }
    
    if (!is.null(markers)){
      ids <- intersect(unlist(markers), rownames(x))
      m.ids = match(ids, rownames(x))
      Sigma <- Sigma[m.ids, ]
    }
    if(verbose){message("Creating Variance Matrix...")}
  }
  
  S <- sapply(unique(clusters), function(ct){
    my.rowMeans(sapply(unique(samples), function(sid){
      y = counts(x)[, clusters %in% ct & samples %in% sid]
      if(is.null(dim(y))){
        return(sum(y))
      }else{
        return(sum(y)/ncol(y))
      }
    }), na.rm = TRUE)
  })
  if(verbose){message("Creating Library Size Matrix...")}
  
  S[S == 0] = NA
  M.S = colMeans(S, na.rm = TRUE)
  #S.ra = relative.ab(S, by.col = FALSE)
  #S.ra[S.ra == 0] = NA
  #S[S == 0] = NA
  #M.S = mean(S, na.rm = TRUE)*ncol(S)*colMeans(S.ra, na.rm = TRUE)
  
  if(!is.null(cell_size)){
    if(!is.data.frame(cell_size)){
      stop("cell_size paramter should be a data.frame with 1st column for cell type names and 2nd column for cell sizes")
    }else if(sum(names(M.S) %in% cell_size[, 1]) != length(names(M.S))){
      stop("Cell type names in cell_size must match clusters")
    }else if (any(is.na(as.numeric(cell_size[, 2])))){
      stop("Cell sizes should all be numeric")
    }
    my_ms_names <- names(M.S)
    cell_size <- cell_size[my_ms_names %in% cell_size[, 1], ]
    M.S <- cell_size[match(my_ms_names, cell_size[, 1]),]
    M.S <- M.S[, 2]
    names(M.S) <- my_ms_names
  }
  
  D <- t(t(M.theta)*M.S)
  
  if(!is.null(select.ct)){
    m.ct = match(select.ct, colnames(D))
    D = D[, m.ct]
    S = S[, m.ct]
    M.S = M.S[m.ct]
    M.theta = M.theta[, m.ct]
  }
  
  if (!is.null(markers)){
    ids <- intersect(unlist(markers), rownames(x))
    m.ids = match(ids, rownames(x))
    D <- D[m.ids, ]
    M.theta <- M.theta[m.ids, ]
  }
  
  if(ct.cov){
    return(list(Disgn.mtx = D, S = S, M.S = M.S, M.theta = M.theta, Sigma.ct = Sigma.ct))
  }else{
    return(list(Disgn.mtx = D, S = S, M.S = M.S, M.theta = M.theta, Sigma = Sigma))
  }
}

music_prop.cluster = function(bulk.mtx, sc.sce, group.markers, groups, clusters, samples, clusters.type,
                              verbose = TRUE, iter.max = 1000, nu = 0.0001, eps = 0.01, centered = FALSE, normalize = FALSE, ... ){
  bulk.gene = rownames(bulk.mtx)[rowMeans(bulk.mtx) != 0]
  bulk.mtx = bulk.mtx[bulk.gene, ]
  select.ct = unlist(clusters.type)
  
  if(length(setdiff(names(group.markers), names(clusters.type))) > 0 || length(setdiff(names(clusters.type), names(group.markers))) > 0){
    stop("Cluster number is not matching!")
  }else{
    group.markers = group.markers[names(clusters.type)]
  }
  
  
  if(verbose){message('Start: cluster estimations...')}
  cluster.sc.basis = music_basis(sc.sce, non.zero = TRUE, markers = NULL, clusters = groups, samples = samples, select.ct = names(clusters.type), verbose = verbose)
  if(verbose){message('Start: cell type estimations...')}
  sc.basis = music_basis(sc.sce, non.zero = TRUE, markers = NULL, clusters = clusters, samples = samples, select.ct = select.ct, verbose = verbose)
  cm.gene = intersect(rownames(sc.basis$Disgn.mtx), bulk.gene)
  
  if(length(cm.gene)< 0.2*min(length(bulk.gene), nrow(sc.sce)) ){
    stop("Too few common genes!")
  }
  
  if(verbose){message(paste('Used', length(cm.gene), 'common genes...'))}
  
  m.sc = match(cm.gene, rownames(sc.basis$Disgn.mtx)); m.bulk = match(cm.gene, bulk.gene)
  group.markers = lapply(group.markers, intersect, cm.gene)
  
  D1 = sc.basis$Disgn.mtx[m.sc,]; M.S = sc.basis$M.S; Sigma = sc.basis$Sigma[m.sc, ]
  
  cluster.select = setdiff(rownames(D1), unique(unlist(group.markers)))
  cluster.diff = unique(unlist(group.markers))
  
  D1.cluster = cluster.sc.basis$Disgn.mtx[cluster.select, ]; M.S.cluster = cluster.sc.basis$M.S;
  Yjg = relative.ab(bulk.mtx[m.bulk, ]); N.bulk = ncol(bulk.mtx);
  
  Sigma.cluster = cluster.sc.basis$Sigma[cluster.select, ];
  
  D1.sub = cluster.sc.basis$Disgn.mtx[cluster.diff, ]; Sigma.sub = cluster.sc.basis$Sigma[cluster.diff, ];
  
  Est.prop.weighted.cluster = NULL
  for(i in 1:N.bulk){
    if(sum(Yjg[, i] == 0) > 0){
      name.temp = rownames(Yjg)[Yjg[, i] != 0]
      D1.cluster.temp = D1.cluster[rownames(D1.cluster)%in%name.temp, ];
      D1.sub.temp = D1.sub[rownames(D1.sub)%in%name.temp, ];
      Yjg.temp = Yjg[Yjg[, i]!=0, i];
      Sigma.cluster.temp = Sigma.cluster[rownames(Sigma.cluster)%in%name.temp, ];
      Sigma.sub.temp = Sigma.sub[rownames(Sigma.sub)%in%name.temp, ];
      if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...') )
    }else{
      D1.cluster.temp = D1.cluster;
      D1.sub.temp = D1.sub;
      Yjg.temp = Yjg[, i];
      Sigma.cluster.temp = Sigma.cluster;
      Sigma.sub.temp = Sigma.sub;
      if(verbose) message(paste(colnames(Yjg)[i], 'has common genes', sum(Yjg[, i] != 0), '...'))
    }
    lm.D1.cluster = music.iter(Yjg.temp, D1.cluster.temp, M.S.cluster, Sigma.cluster.temp, iter.max = iter.max,
                               nu = nu, eps = eps, centered = centered, normalize = normalize)
    p.weight = NULL
    p.cluster.weight = lm.D1.cluster$p.weight
    for(j in 1:length(clusters.type)){
      if(length(clusters.type[[j]]) == 1){
        p.weight = c(p.weight, p.cluster.weight[j])
      }else{
        if(p.cluster.weight[j] == 0){
          p.weight = c(p.weight, rep(0, length(clusters.type[[j]])))
        }else{
          c.marker = intersect(group.markers[[j]], names(Yjg.temp))
          Y.sub = D1.sub.temp[c.marker, j]*p.cluster.weight[j] +
            (Yjg.temp[c.marker] - D1.sub.temp[c.marker, ]%*% p.cluster.weight) * p.cluster.weight[j]
          names(Y.sub) = c.marker
          Y.sub = Y.sub[Y.sub > 0]
          lm.D1.sub = music.iter(Y.sub, D1[c.marker, clusters.type[[j]]], M.S[clusters.type[[j]]],
                                 Sigma[c.marker, clusters.type[[j]]])
          p.weight = c(p.weight, p.cluster.weight[j] * lm.D1.sub$p.weight)
        }
      }
    }
    Est.prop.weighted.cluster = rbind(Est.prop.weighted.cluster, p.weight)
  }
  
  colnames(Est.prop.weighted.cluster) = unlist(clusters.type)
  rownames(Est.prop.weighted.cluster) = colnames(Yjg)
  
  return(list(Est.prop.weighted.cluster = Est.prop.weighted.cluster))
}



