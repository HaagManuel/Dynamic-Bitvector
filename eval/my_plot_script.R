# Loads libraries
library(tidyverse)
library(dplyr)
library(scales)
library(ggpubr)
library(ggtext)
library(scales)


setwd("C:/Users/manue/OneDrive/Dokumente/uni/master/2.Semester/Fortgeschrittene Datenstrukturen/project/experiments")
out_folder='./'
header1 = c("algo", "run", "n", "leaf_blocks", "insert", "remove", "memory_leafs", "memory_nodes", "memory_total", "fraction_stored_used_bits","bits_per_stored_bit") 
header2 = c("algo", "run", "n", "insert") 
header3 = c("algo", "run", "n", "leaf_blocks", "flip")
header4 = c("algo", "run", "n", "leaf_blocks", "front", "back", "rec", "compress_rec")

add_group_label = function(df) {
  df$Group = paste(df$algo)
  return(df)
}

build_up_big = read.table("./build_up_bv_big.csv", col.names = header1, comment.char = '#')
build_up_big_oop = read.table("./build_up_bv_oop_big.csv", col.names = header2, comment.char = '#')
build_up_big_oop$leaf_blocks = 64
build_up_big_bp = read.table("./build_up_bp_big.csv", col.names = header1, comment.char = '#')

build_up = read.table("./build_up_bv.csv", col.names = header1, comment.char = '#')
build_up_oop = read.table("./build_up_bv_oop.csv", col.names = header2, comment.char = '#')
by_leafs_small = read.table("./by_leafs.csv", col.names = header1, comment.char = '#')
by_leafs_big = read.table("./by_leafs_big.csv", col.names = header1, comment.char = '#')

construction_bv = read.table("./construction_bv.csv", col.names = header4, comment.char = '#')

by_leafs = by_leafs_big

flip_bv = read.table("./flip_bv.csv", col.names = header3, comment.char = '#')

eps = 0.001

prepare_df <- function(df, time_col, algo_name, select_cols) { 
  new_df = df
  new_df$time = time_col
  new_df$algo = algo_name
  new_df = select(new_df, select_cols)
  new_df
}

time_cols = c("algo", "leaf_blocks", "run", "n", "time")


build_up_insert = prepare_df(build_up_big, build_up_big$insert, "bv_insert", time_cols)
build_up_remove = prepare_df(build_up_big, build_up_big$remove, "bv_remove", time_cols)
build_up_insert_oop = prepare_df(build_up_big_oop, build_up_big_oop$insert, "bv_insert_oop", time_cols)
build_up_insert_bp = prepare_df(build_up_big_bp, build_up_big_bp$insert, "bp_insert", time_cols)
build_up_remove_bp = prepare_df(build_up_big_bp, build_up_big_bp$remove, "bp_remove", time_cols)


exp1_df = rbind(build_up_insert, build_up_remove, build_up_insert_oop, build_up_insert_bp,build_up_remove_bp )
exp1_df = add_group_label(exp1_df)
exp1_df[, "time"] = exp1_df[, "time"] + eps #avoid 0 in log


leaf_insert = by_leafs
leaf_insert$time = by_leafs$insert
leaf_insert$algo = paste(by_leafs$algo, "insert", sep="_")
leaf_insert = select(leaf_insert, c("algo", "run", "leaf_blocks", "time"))

leaf_remove = by_leafs
leaf_remove$time = leaf_remove$remove
leaf_remove$algo = paste(by_leafs$algo, "remove", sep="_")
leaf_remove = select(leaf_remove, c("algo", "run", "leaf_blocks", "time"))

exp_leaf = rbind(leaf_insert, leaf_remove)
exp_leaf = add_group_label(exp_leaf)

exp_memory_leaf = by_leafs
exp_memory_leaf = add_group_label(exp_memory_leaf)

exp_flip = flip_bv
exp_flip$time = exp_flip$flip
exp_flip = add_group_label(exp_flip)



construction_front =  prepare_df(construction_bv, construction_bv$front, "front", time_cols)
construction_back =  prepare_df(construction_bv, construction_bv$back, "back", time_cols)
construction_rec =  prepare_df(construction_bv, construction_bv$rec, "rec", time_cols)
construction_compr_rec =  prepare_df(construction_bv, construction_bv$compress_rec, "compress_rec", time_cols)
exp_construction = rbind(construction_front, construction_back, construction_rec, construction_compr_rec)
exp_construction = add_group_label(exp_construction)
exp_construction[, "time"] = exp_construction[, "time"] + eps #avoid 0 in log

############## Plot time n ####################
# This plots the mean and standard error for each section, and connects the points with lines
plot_time_n <- function(data, title) {
  num_variants = length(unique(data$Group))
  max_k = log2(max(data$n))
  #ns = 2^(1:max_k)
  ns = 2^(2*(1:(max_k/2)))
  labels = paste("2^", 1:max_k, sep="")
  
  min_k_time = floor(log2(min(data$time)))
  max_k_time = ceiling(log2(max(data$time)))
  #ns_time = 2^(min_k_time:max_k_time)
  ns_time = 2^(2*((min_k_time/2):(max_k/2)))
  labels_time = paste("2^", min_k_time:max_k_time, sep="")

  ggplot(data, aes(x = n, y = time, group = Group, color = Group, shape = Group)) +
    stat_summary(fun.data = mean_se, geom = 'pointrange') +
    stat_summary(fun = mean, geom = 'line') +
    scale_x_continuous(trans="log2", breaks=ns, labels=trans_format("log2", math_format(2^.x))) +
    scale_y_continuous(trans="log2", breaks=ns_time, labels=trans_format("log2", math_format(2^.x))) +
    scale_shape_manual(values=1:num_variants) +
    ggtitle(title) +
    theme(text = element_text(size = 20)) +
    xlab("n") +
    ylab("Avg Time [ms]")
    
}

plot_time_n(exp1_df, "Average runtime by number of operations n")
plot_time_n(exp_flip, "Average runtime of 10e6 flips at every n ")
plot_time_n(subset(exp_construction, leaf_blocks == 64), "Average construction time")
############## Plot time n ####################


############## Plot memory ####################
plot_memory <- function(data, title="") {
  # This plots the mean and standard error for each section, and connects the points with lines
  num_variants = length(unique(data$Group))
  ggplot(data, aes(x = leaf_blocks, y = bits_per_stored_bit, group = Group, color = Group, shape = Group)) +
    ylim(1, 4) +
    stat_summary(fun.data = mean_se, geom = 'pointrange') +
    stat_summary(fun = mean, geom = 'line') +
    scale_x_continuous(trans="log2", breaks=unique(data$leaf_blocks)) +
    scale_y_continuous(breaks=seq(1,4, 0.5)) +
    scale_shape_manual(values=1:num_variants) +
    ggtitle(title) +
    xlab("Leaf Blocks") +
    ylab("Bits per stored bit") +
    theme(text = element_text(size = 20)) +
    geom_hline(yintercept=1)
}

plot_memory(exp_memory_leaf, "")

############## Plot memory ####################

############## Plot time by leaf ####################
plot_time_leaf <- function(data, title="") {
  # This plots the mean and standard error for each section, and connects the points with lines
  min_k_time = floor(log2(min(data$time)))
  max_k_time = ceiling(log2(max(data$time)))
  #ns_time = 2^(min_k_time:max_k_time)
  ns_time = 2^(2*((min_k_time/2):(max_k_time/2)))
  labels_time = paste("2^", min_k_time:max_k_time, sep="")
  
  num_variants = length(unique(data$Group))
  ggplot(data, aes(x = leaf_blocks, y = time, group = Group, color = Group, shape = Group)) +
    stat_summary(fun.data = mean_se, geom = 'pointrange') +
    stat_summary(fun = mean, geom = 'line') +
    scale_x_continuous(trans="log2", breaks=unique(data$leaf_blocks)) +
    #scale_y_continuous(trans="log2", breaks=ns_time, labels=labels_time) +
    scale_y_continuous(trans="log2", breaks=ns_time, labels=trans_format("log2", math_format(2^.x))) +
    scale_shape_manual(values=1:num_variants) +
    ggtitle(title) +
    theme(text = element_text(size = 20)) +
    xlab("Leaf Blocks") +
    ylab("Avg Time [ms]")
}
#plot_time_leaf(exp_leaf, "Average runtime by leaf sizes, n=2^18")
plot_time_leaf(exp_leaf, bquote("Average runtime by leaf sizes, n =" ~ 2^18))
plot_time_leaf(subset(exp_construction, n == max(exp_construction$n)), "Construction time by leaf sizes")


#bquote("Hello" ~ r[xy] == .(2) ~ "and" ~ B^2)
############## Plot time by leaf ####################





#run time insert, remove
p1 = plot_time_n(exp1_df, "run time by number of operations n")

#flip operation
p2 = plot_time_n(exp_flip, bquote("run time of m =" ~ 10^6 ~ "flip operations"))

#run time by leafs
p3 = plot_time_leaf(exp_leaf, bquote("run time by leaf sizes, n =" ~ 2^18))

#memory by leafs
p4 = plot_memory(exp_memory_leaf, bquote("space usage by leaf sizes, n =" ~ 2^18))

#construction
p5 = plot_time_n(subset(exp_construction, leaf_blocks == 64), "construction time, 64 leaf blocks")
p6 = plot_time_leaf(subset(exp_construction, n == max(exp_construction$n)), "construction time by leaf sizes")

p1
p2
p3
p4
p5
p6


plot_height = 467 
plot_width = 788


############## PNG  ####################
png(paste(out_folder, "insert_remove_by_n.png", sep = "/"), width = plot_width, height = plot_height)
p1
dev.off()

png(paste(out_folder, "flip_by_n.png", sep = "/"), width = plot_width, height = plot_height)
p2
dev.off()

png(paste(out_folder, "time_by_leaf.png", sep = "/"), width = plot_width, height = plot_height)
p3
dev.off()

png(paste(out_folder, "memory_by_leaf.png", sep = "/"), width = plot_width, height = plot_height)
p4
dev.off()

png(paste(out_folder, "construction_by_n.png", sep = "/"), width = plot_width, height = plot_height)
p5
dev.off()

png(paste(out_folder, "construction_by_leaf.png", sep = "/"), width = plot_width, height = plot_height)
p6
dev.off()

