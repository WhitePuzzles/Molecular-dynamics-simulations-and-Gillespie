#能量最小化 采用梯度下降
export min_energy



function min_energy(step_tol,F_tol,part_num,initial_pos,mass,force_types,boundary)
    λ=0.1
    force=zeros(part_num,3)
    force_sca=zeros(part_num,1)
    pos=copy(initial_pos)
    step=0
    f_max=Inf
    #f_temp=zeros(3)
    #dist=zeros(3)
    thread_num=Threads.nthreads() #多线程
    force_chunk=[zero(force) for i in 1:thread_num]
    while step<step_tol && f_max>F_tol
        force.=0.0
        force_sum!(force,part_num,force_types,mass,pos, force_chunk,thread_num)
        force_sca=sqrt.(sum(force.^2,dims=2))
        pos.=pos.+λ*force./force_sca
        step+=1
        f_max=maximum(force_sca)
    end
    return pos
end