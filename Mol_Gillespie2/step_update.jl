#运动迭代函数
export step_update!


#运动迭代函数
function step_update!(part_num,mass,pos,vel,force_types,force,boundary,f_temp,iterate_method,other_para)
    if iterate_method=="velocity verlet" #Velocity Verlet integration method
        @. pos=pos+vel+force/mass/2
        @. vel=vel+force/mass/2
        force.=0.0
        for force_i in force_types
            add_force!(part_num,mass,pos,force_i,force,f_temp)
        end
        @. vel=vel+force/mass/2
    elseif iterate_method=="leap-frog verlet" #leap-frog Verlet integration method
        @. vel=vel+force/mass
        @. pos=pos+vel
        force.=0.0
        for force_i in force_types
            add_force!(part_num,mass,pos,force_i,force,f_temp)
        end
    elseif iterate_method=="langevin leap-frog" #langevin leap-frog integration method
        @. pos=pos+vel
        force.=0.0
        for force_i in force_types
            add_force!(part_num,mass,pos,force_i,force,f_temp)
        end
        f_temp2=randn(part_num,3)
        @. vel=vel*other_para+(1-other_para)/γ*force/mass+sqrt(kB*T*(1-other_para^2))*f_temp2/sqrt(mass)
    end
end