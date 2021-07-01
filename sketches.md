poissonDisc_sample(分散均匀)



Type_b;

Type_p;

浓度（基底的粒子c_b

c_p : 粒子浓度

attract force (type...; 可以注册作用力(1/r2 1/r3 1/r，$\alpha$，或者其他作用力, 作用范围)

结合概率 p

结合范围 r_e

> 1:判断是否结合上: 结合有1-p概率逃脱，单次更新的步长不能太长
>
> 2:粒子之间有体积力(防止结合)
>
> 3:随机取向($\sigma$)

更新过程 通过id索引 计算范围R(矩形范围) 

通过排序 避免O(N^2)的复杂度，k N logN；









