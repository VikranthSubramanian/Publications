clc;
clear all;
close all;

f1=ode23('lpnnnnn',[0 0.005],[100 200 rand(3,1)' rand(3,1)' rand(3,1)' rand(3,1)' 0 0 0 0 0 0 0 0 0 0 0 0 0 0]');

