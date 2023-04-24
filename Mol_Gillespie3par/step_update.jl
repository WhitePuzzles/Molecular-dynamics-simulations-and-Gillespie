#运动迭代函数
export step_update!



function step_update!(part_num,mass,pos,vel,force_types,force,boundary,iterate_method,other_para,force_chunk,thread_num)
    if iterate_method=="velocity verlet" #Velocity Verlet integration method
        @. pos=pos+vel+force/mass/2
        @. vel=vel+force/mass/2
        force_sum!(force,part_num,force_types,mass,pos,force_chunk,thread_num)
        @. vel=vel+force/mass/2
    elseif iterate_method=="leap-frog verlet" #leap-frog Verlet integration method
        @. vel=vel+force/mass
        @. pos=pos+vel
        force_sum!(force,part_num,force_types,mass,pos,force_chunk,thread_num)
    elseif iterate_method=="langevin leap-frog" #langevin leap-frog integration method
        @. pos=pos+vel
        force_sum!(force,part_num,force_types,mass,pos,force_chunk,thread_num)
        rand_n=randn(part_num,3)
        @. vel=vel*other_para+(1-other_para)/γ*force/mass+sqrt(kB*T*(1-other_para^2))*rand_n/sqrt(mass)
    end
end