#力场
export
    Lennard_Jones,
    WCA,
    Rouse_model,
    FENE,
    Coulomb,
    add_force,
    force_sum!




#累加函数
function force_sum!(force,part_num,force_types,mass,pos,force_chunk,thread_num)
    for i in 1:thread_num
        force_chunk[i].=0.
    end
    Threads.@threads for thread_id in 1:thread_num
        for p_i in thread_id:thread_num:(part_num-1)
            for p_j in (p_i+1):part_num
                for force_type in force_types
                    f_temp=add_force(mass,pos,force_type,p_i,p_j)
                    force_chunk[thread_id][p_i,:]+=f_temp
                    force_chunk[thread_id][p_j,:]-=f_temp
                end
            end
        end
    end
    force.=sum(force_chunk)
end



#力场
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
function add_force(mass,pos,force_type::Lennard_Jones,p_i,p_j)
    return -24*force_type.epsilon*force_type.sigma^6/norm(pos[p_j,:]-pos[p_i,:])^8*(2*force_type.sigma^6/norm(pos[p_j,:]-pos[p_i,:])^6-1)*(pos[p_j,:]-pos[p_i,:])
end

#WCA
function add_force(mass,pos,force_type::WCA,p_i,p_j)
    return norm(pos[p_j,:]-pos[p_i,:])<2^(1/6)*force_type.sigma ? -24*force_type.epsilon*force_type.sigma^6/norm(pos[p_j,:]-pos[p_i,:])^8*(2*force_type.sigma^6/norm(pos[p_j,:]-pos[p_i,:])^6-1)*(pos[p_j,:]-pos[p_i,:]) : [0.,0.,0.]
end

#Rouse_model
function add_force(mass,pos,force_type::Rouse_model,p_i,p_j)
    return force_type.M[p_i,p_j]>0 ? force_type.M[p_i,p_j]*(1-force_type.L0[p_i,p_j]/norm(pos[p_j,:]-pos[p_i,:]))*(pos[p_j,:]-pos[p_i,:]) : [0.,0.,0.]
end

#finite extensible nonlinear elastic(FENE)
function add_force(mass,pos,force_type::FENE,p_i,p_j)
    return (force_type.M[p_i,p_j]>0 && norm(dist)<force_type.R0) ? force_type.M[p_i,p_j]/(force_type.sigma^2(1-(norm(pos[p_j,:]-pos[p_i,:])/force_type.R0)^2))*(pos[p_j,:]-pos[p_i,:]) : [0.,0.,0.]
end

#Coulomb
function add_force(mass,pos,force_type::Coulomb,p_i,p_j)
    return -kC*force_type.V[p_i]*force_type.V[p_j]/norm(pos[p_j,:]-pos[p_i,:])^2*normalize(pos[p_j,:]-pos[p_i,:])
end