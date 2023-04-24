#力场
export
    Lennard_Jones,
    WCA,
    Rouse_model,
    FENE,
    Coulomb,
    add_force!



struct Lennard_Jones #U=4ε((σ/r)^12-(σ/r)^6)
    epsilon::Float64
    sigma::Float64
end

struct WCA #U=4ε((σ/r)^12-(σ/r)^6+1/4),r<2^{1/6}σ. 0,r>2^{1/6}σ
    epsilon::Float64
    sigma::Float64
end

struct Rouse_model #每一对的能量U=1/2*k*(r-L0)^2. M是邻接矩阵，ij粒子连接且弹性系数为k_ij=M[i,j]，M[i,j]=-1表示不连接.L0是相邻两粒子平均距离
    M::Matrix
    L0::Matrix
end

struct FENE #U=-1/2k*R0^2/σ^2*ln(1-(r/R0)^2),r≤R0. 0,r>R0
    M::Matrix
    R0::Float64
    sigma::Float64
end

struct Coulomb #U=kq1q2/r
    V::Vector #storage charge quantity
end

#Lennard_Jones
function add_force!(part_num,mass,pos,force_type::Lennard_Jones,force,f_temp)
    pos_dist=zeros(part_num,3)
    for i in 1:part_num-1
        pos_dist.=pos.-transpose(pos[i,:])
        for j in (i+1):part_num
            f_temp.=-24*force_type.epsilon*force_type.sigma^6/norm(@view pos_dist[j,:])^8*(2*force_type.sigma^6/norm(@view pos_dist[j,:])^6-1)* (@view pos_dist[j,:])
            force[i,:]+=f_temp
            force[j,:]-=f_temp
        end
    end
end

#WCA
function add_force!(part_num,mass,pos,force_type::WCA,force,f_temp)
    pos_dist=zeros(part_num,3)
    for i in 1:part_num-1
        pos_dist.=pos.-transpose(pos[i,:])
        for j in (i+1):part_num
            f_temp.=norm(@view pos_dist[j,:])<2^(1/6)*force_type.sigma ? -24*force_type.epsilon*force_type.sigma^6/norm(@view pos_dist[j,:])^8*(2*force_type.sigma^6/norm(@view pos_dist[j,:])^6-1)* (@view pos_dist[j,:]) : [0.,0.,0.]
            force[i,:]+=f_temp
            force[j,:]-=f_temp
        end
    end
end

#Rouse_model
function add_force!(part_num,mass,pos,force_type::Rouse_model,force,f_temp)
    pos_dist=zeros(part_num,3)
    for i in 1:part_num-1
        pos_dist.=pos.-transpose(pos[i,:])
        for j in (i+1):part_num
            if force_type.M[i,j]>0
                #f_temp.=force_type.M[i,j]*(1-force_type.L0[i,j]/norm(pos[j,:]-pos[i,:]))*(pos[j,:]-pos[i,:])
                f_temp.=force_type.M[i,j]*(1-force_type.L0[i,j]/norm(@view pos_dist[j,:]))*(@view pos_dist[j,:])
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
        pos_dist.=pos.-transpose(pos[i,:])
        for j in (i+1):part_num
            if force_type.M[i,j]>0
                #f_temp.= norm(pos[j,:]-pos[i,:])<force_type.R0 ? force_type.M[i,j]/(force_type.σ^2(1-(norm(pos[j,:]-pos[i,:])/force_type.R0)^2))*(pos[j,:]-pos[i,:]) : [0,0,0]
                f_temp.= norm(@view pos_dist[j,:])<force_type.R0 ? force_type.M[i,j]/(force_type.sigma^2(1-(norm(@view pos_dist[j,:])/force_type.R0)^2))*(@view pos_dist[j,:]) : [0.,0.,0.]
                force[i,:]+=f_temp
                force[j,:]-=f_temp
            end
        end
    end
end

function add_force!(part_num,mass,pos,force_type::Coulomb,force,f_temp)
    pos_dist=zeros(part_num,3)
    for i in 1:part_num-1
        pos_dist.=pos.-transpose(pos[i,:])
        for j in (i+1):part_num
            f_temp.=-kC*force_type.V[i]*force_type.V[j]/norm(@view pos_dist[j,:])^2*normalize(@view pos_dist[j,:])
            force[i,:]+=f_temp
            force[j,:]-=f_temp
        end
    end
end