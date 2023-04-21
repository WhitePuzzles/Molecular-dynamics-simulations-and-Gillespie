#力场
export
    Lennard_Jones,
    Rouse_model,
    FENE,
    add_force!





struct Lennard_Jones #U=4ε((σ/r)^12-(σ/r)^6)
    ε::Float64
    σ::Float64
end

struct Rouse_model #每一对的能量U=1/2*k*(r-L0)^2. M是邻接矩阵，ij粒子连接且弹性系数为k_ij=M[i,j]，M[i,j]=-1表示不连接.L0是相邻两粒子平均距离
    M::Matrix
    L0::Matrix
end

struct FENE #U=-1/2k*R0^2/σ^2*ln(1-(r/R0)^2),r≤R0. 0,r>R0
    M::Matrix
    R0::Float64
    σ::Float64
end
#Lennard_Jones
function add_force!(part_num,mass,pos,force_type::Lennard_Jones,force,f_temp)
    pos_dist=zeros(part_num,3)
    for i in 1:part_num-1
        @. pos_dist=pos-pos[i,:]'
        for j in (i+1):part_num
            f_temp.=-24*force_type.ε*force_type.σ^6/norm(pos_dist[j,:])^8*(2*force_type.σ^6/norm(pos_dist[j,:])^6-1)*pos_dist[j,:]
            force[i,:]+=f_temp
            force[j,:]-=f_temp
        end
    end
end

#Rouse_model
function add_force!(part_num,mass,pos,force_type::Rouse_model,force,f_temp)
    pos_dist=zeros(part_num,3)
    for i in 1:part_num-1
        for j in (i+1):part_num
            if force_type.M[i,j]>0
                f_temp.=force_type.M[i,j]*(1-force_type.L0[i,j]/norm(pos[j,:]-pos[i,:]))*(pos[j,:]-pos[i,:])
                force[i,:]+=f_temp
                force[j,:]-=f_temp
            end
        end
    end
end

#finite extensible nonlinear elastic(FENE)
function add_force!(part_num,mass,pos,force_type::FENE,force,f_temp)
    pos_dist=zeros(part_num,3)
    for i in 1:part_num-1
        for j in (i+1):part_num
            if force_type.M[i,j]>0
                f_temp.= norm(pos[j,:]-pos[i,:])<force_type.R0 ? force_type.M[i,j]/(force_type.σ^2(1-(norm(pos[j,:]-pos[i,:])/force_type.R0)^2))*(pos[j,:]-pos[i,:]) : [0,0,0]
                force[i,:]+=f_temp
                force[j,:]-=f_temp
            end
        end
    end
end