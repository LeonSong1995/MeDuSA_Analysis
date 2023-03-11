#! /usr/lib/R/bin/Rscript
library(tidyr)
library(argparser, quietly=TRUE)

p = arg_parser("simulation")
p = add_argument(p, "prefix", help="")
p = add_argument(p, "resolution", help="")
p = add_argument(p, "dataPath", help="")
p = add_argument(p, "OutPath", help="")
argv = parse_args(p)


# Generate pseuod-bulk using the selected cells
generate_pseudo_bulk = function(abundance,val_mat,cellLocationOnGrid){
  num = round(1/min(abundance[which(abundance>0)]))
  used_state = unique(cellLocationOnGrid)
  used_state = used_state[order(readr::parse_number(used_state))]

  #sample cells based on the abundance of each cell-state
  cell = unlist(lapply(used_state, function(state) {
    print(state)
    index = names(which(cellLocationOnGrid == state))
    sample(index, round(num * abundance[match(state, used_state)]), replace = TRUE)
  }))
    
  #average the expression profiles of selected cells 
  chunk_size = 10000
  n_chunks = ceiling(length(cell) / chunk_size)

  pseudo_bulk = numeric(nrow(val_mat))
  for (i in seq_len(n_chunks)) {
    print(paste0('merging cell of chunk:',i))
    start =(i - 1) * chunk_size + 1
    end = min(i * chunk_size, length(cell))
    subset = cell[start:end]
    pseudo_bulk= pseudo_bulk + rowSums(val_mat[,subset])
  }
  pseudo_bulk = pseudo_bulk / length(cell)


  return(pseudo_bulk)
}

# Simulation function (Select Cells)
simulation_cell = function(traj_ref,traj_val,val_mat,k=4,resolution=50){

  #binning the trajectory
  traj_ref = as.matrix(traj_ref / max(traj_ref))
  traj_val = as.matrix(traj_val / max(traj_val))

  #sort the traj_val by increasing
  traj_val = traj_val[order(traj_val,decreasing=F),,drop=F]
  val_mat = val_mat[,rownames(traj_val)]

  #cut the cell-state trajectory in to bins
  bin =  paste0('bin',cut(traj_ref, resolution, labels = FALSE, include.lowest = TRUE))
  breaks = aggregate(traj_ref,by=list(bin),FUN=min)
  breaks = breaks[order(readr::parse_number(breaks[,1])),2]
  breaks[1]=-Inf
  breaks[length(breaks)+1]=Inf


  #assin cells in simulation sets to each cell state. 
  cellLocationOnGrid = rep(NA,length(traj_val))
  for(currBreakIndex in 2:length(breaks)){cellLocationOnGrid[which(traj_val >=breaks[currBreakIndex-1] & traj_val <breaks[currBreakIndex])]=currBreakIndex-1}
  cellLocationOnGrid = paste0('bin',cellLocationOnGrid)
  names(cellLocationOnGrid) = rownames(traj_val)


  #median pseudo-time of each cell-state bin 
  bmed = aggregate(traj_ref,by=list(bin),FUN=median)
  rownames(bmed) = bmed[,1]
  bmed = bmed[order(readr::parse_number(bmed[,1])),2]


  #generate pseudo bulk RNA-seq data 
  #shape-1: Decreasing
  ab1 =  -(bmed)^1+abs(min(-(bmed)^1))
  ab1 = ab1^k
  ab1 = ab1/sum(ab1)
  b1 = generate_pseudo_bulk(ab1,val_mat,cellLocationOnGrid)
  print('Decreasing')

  #shape-2: Incremental
  ab2 = (bmed)^1
  ab2 = ab2^k
  ab2 = ab2/sum(ab2)
  b2 = generate_pseudo_bulk(ab2,val_mat,cellLocationOnGrid)
  print('Increasing')

  #shape-3: Unimodal
  ab3 = -((bmed-bmed[round(length(bmed)/2)])^2)
  ab3 = ab3 + max(((bmed-bmed[round(length(bmed)/2)])^2))
  ab3 = ab3^k
  ab3 =  ab3/sum(ab3)
  b3 = generate_pseudo_bulk(ab3,val_mat,cellLocationOnGrid)
  print('Unimodal')

  #shape-4: Bimodal
  ab4 = sin(3*pi*(bmed/max(bmed)))
  ab4 = ab4 -min(ab4)
  ab4= ab4^1
  ab4 =  ab4/sum(ab4)
  b4 = generate_pseudo_bulk(ab4,val_mat,cellLocationOnGrid)
  print('Bimodal')

  bulk = cbind(b1,b2,b3,b4)
  abundace = cbind(ab1,ab2,ab3,ab4)
  colnames(bulk) = colnames(abundace) = c('Decreasing','Increasing','Unimodal','Bimodal')
  rownames(abundace) = paste0('bin',seq(1,nrow(abundace)))

  return(list(bulk,abundace))
}

# Simulation function (Mean Cells)
simulation_mean = function(traj_ref,traj_val,val_mat,k=4,resolution=50){

  #binning the trajectory
  traj_ref = as.matrix(traj_ref / max(traj_ref))
  traj_val = as.matrix(traj_val / max(traj_val))

  bin = MeDuSA::Partition_cell_trajectory(traj_ref,nbins = resolution)
  breaks = aggregate(traj_ref,by=list(bin),FUN=min)
  breaks = breaks[order(readr::parse_number(breaks[,1])),2]
  breaks[1]=-Inf
  breaks[length(breaks)+1]=Inf

  # generate pseudo bulk RNA-seq data using other samples####
  traj_val = sort(traj_val[,1]) / max(traj_val)
  val_mat = val_mat[,names(traj_val)]


  #assin cells in simulation sets to each cell state. 
  cellLocationOnGrid = rep(NA,length(traj_val))
  for(currBreakIndex in 2:length(breaks)){cellLocationOnGrid[which(traj_val >=breaks[currBreakIndex-1] & traj_val <breaks[currBreakIndex])]=currBreakIndex-1}
  names(cellLocationOnGrid)=colnames(val_mat)


  #generate the cell-state gene expression matrix 
  refz = sapply(unique(cellLocationOnGrid),function(j){
    index = which(cellLocationOnGrid==j)
    if(length(index)>1){
      out = rowMeans(val_mat[,index])
    }else if (length(index)==1){
      out = val_mat[,index]
    }else{
      print('NOTION!')
      out = 0
    }
    return(out)
  })

  #median pseudo-time of each cell-state bin 
  used_state = unique(bin)
  used_state = used_state[order(readr::parse_number(used_state))]
  bmed = aggregate(traj_ref,by=list(bin),FUN=median)
  rownames(bmed) = bmed[,1]
  bmed = bmed[used_state,2]

  #generate the pseudo-bulk dara
  #shape1: Decreasing
  ab1 =  -(bmed)^1+abs(min(-(bmed)^1))
  ab1 = ab1^k
  ab1 = ab1/sum(ab1)
  b1 =(refz %*% ab1)

  #shape2: Increasing
  ab2 = (bmed)^1
  ab2 = ab2^k
  ab2 =  ab2/sum(ab2)
  b2 =(refz %*% ab2)

  #shape3: Unimodal
  ab3 = -((bmed-bmed[round(length(bmed)/2)])^2)
  ab3 = ab3 + max(((bmed-bmed[round(length(bmed)/2)])^2))
  ab3 = ab3^k
  ab3 =  ab3/sum(ab3)
  b3 =(refz %*% ab3)

  #shape4: Bimodal
  ab4 = sin(3*pi*(bmed/max(bmed)))
  ab4 = ab4 -min(ab4)
  ab4= ab4^1
  ab4 =  ab4/sum(ab4)
  b4 =(refz %*% ab4)

  bulk = cbind(b1,b2,b3,b4)
  abundance = cbind(ab1,ab2,ab3,ab4)
  colnames(bulk) = colnames(abundance) = c('Decreasing','Increasing','Unimodal','Bimodal')
  rownames(abundance) = used_state
  return(list(bulk,abundance,bmed))
}


# Load the data
data = readRDS(paste0(argv$dataPath,argv$prefix))
score = data$cell_trajectory
count = as.matrix(data@assays$RNA@counts)

#Seperate the scRNA-seq data into simulation and validation set
ref_cell = sample(names(score),round(length(score)*0.5))
val_cell = names(score)[!names(score) %in% ref_cell]

traj_ref = score[ref_cell]
traj_val = score[val_cell]

ref_mat = count[,ref_cell]
val_mat = count[,val_cell]

#Run simulation
resolution = as.numeric(argv$resolution)

# Simu = simulation_mean(traj_ref=traj_ref,traj_val=traj_val,val_mat=val_mat,resolution = resolution)
Simu = simulation(traj_ref=traj_ref,traj_val=traj_val,val_mat=val_mat,resolution = resolution)
bulk_simu = Simu[[1]];rownames(bulk_simu) = rownames(count)
abundance_simu = Simu[[2]]

# Add additional noise to the bulk data   
residual = replicate(ncol(bulk_simu),rlnorm(nrow(bulk_simu), 0, noise))
# Remove the extreme value 
residual[which(residual>5*median(residual))] = median(residual)
bulk_simu_noise = bulk_simu + residual
rownames(bulk_simu_noise) = rownames(bulk_simu)

saveRDS(bulk_simu_noise,'pseudo-bulk.rds')
saveRDS(abundance_simu,'simulated-abundance.rds')





