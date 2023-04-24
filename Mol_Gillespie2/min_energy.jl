#能量最小化 采用梯度下降
export min_energy



function min_energy(step_tol,F_tol,part_num,initial_pos,mass,force_types,boundary)
    λ=0.1
    force=zeros(part_num,3)
    f_temp=zeros(3)
    force_sca=zeros(part_num,1)
    pos=copy(initial_pos)
    step=0
    f_max=Inf
    while step<step_tol && f_max>F_tol
        force.=0.0
        for force_i in force_types
            add_force!(part_num,mass,pos,force_i,force,f_temp)
        end
        force_sca=sqrt.(sum(force.^2,dims=2))
        pos.=pos.+λ*force./force_sca
        
        step+=1
        f_max=maximum(force_sca)
    end
    return pos
end