#分子动力学模拟和化学反应耦合
export molecular_simulation_and_Gillespie


#主函数
function molecular_simulation_and_Gillespie(initial_state,react_left,react_right,react_rate,EP_index,distance_affect_rate,  t_max,part_num,mass,initial_pos,initial_vel,force_types,iterate_method,boundary,frame_interval)
    f_mu=zeros(Float64,size(react_left,1))
    f_mu_cumsum=zeros(Float64,size(react_left,1))
    state=copy(initial_state)
    state_record=[copy(state)] #状态记录
    t_record=[0.] #时间节点记录
    rea_tr=transpose(react_left)
    
    t=0
    vel=copy(initial_vel)
    pos=copy(initial_pos)
    force=zeros(part_num,3)
    other_para=0
    frame_record=zeros(part_num,3,1+ceil(Int64,t_max/frame_interval)) #用于收集帧
    count_frame=0
    
    thread_num=Threads.nthreads() #多线程
    force_chunk=[zero(force) for i in 1:thread_num]
    
    #计算初始力
    force_sum!(force,part_num,force_types,mass,pos, force_chunk,thread_num)
    
    #预处理
    if iterate_method=="leap-frog verlet" #预估蛙跳verlet算法初速度
        @. vel=vel-force/mass/2
    elseif iterate_method=="langevin leap-frog" #计算langevin leap-frog算法迭代起始速度
        rand_n0=randn(part_num,3)
        @. vel=vel*exp(-γ/2)+(1-exp(-γ/2))/γ*force/mass+sqrt(kB*T*(1-exp(-γ)))*rand_n0/sqrt(mass)
        alpha=exp(-γ)
        other_para=alpha
    end
    
    
    
    #主循环
    while t<t_max 
        #采样
        if t%frame_interval==0
            count_frame+=1
            frame_record[:,:,count_frame]=pos
        end
            
        #分子动力学
        #更新位置，速度，力
        step_update!(part_num,mass,pos,vel,force_types,force,boundary,iterate_method,other_para, force_chunk,thread_num)
            
        #化学反应
        if t>=t_record[end]
            #f_mu=反应常数+距离影响
            f_mu.=(react_rate+distance_affect_rate(norm(pos[EP_index[1],:]-pos[EP_index[2],:]))).*vec(prod(binomial.(state,rea_tr);dims=1))
            F=sum(f_mu)
            r1=rand()
            tau=-log(r1)/F #产生反应时间
            r2=rand()
            f_mu_cumsum.=cumsum(f_mu)
            next_mu=findfirst(x -> x>F*r2,f_mu_cumsum) #产生反应通道
            state.=state.-react_left[next_mu,:].+react_right[next_mu,:] #状态变化
            push!(state_record,copy(state)) #记录状态
            push!(t_record,t_record[end]+tau) #记录时间节点
        end
        t+=1
    end
    return (pos,state,frame_record,state_record,t_record)
end