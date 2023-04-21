#Gillespie算法
export Gillespie


function Gillespie(stoichiometry_matrix_left,stoichiometry_matrix_right,reaction_rate,initial_state,t_max)
    f_mu=zeros(Float64,size(stoichiometry_matrix_left,1)) #反应倾向函数
    state=initial_state[:]
    state_record=[state[:]] #状态记录
    t=[0.] #时间节点记录
    sml_t=stoichiometry_matrix_left'
    f_mu_cumsum=zeros(Float64,size(stoichiometry_matrix_left,1))
    while t[end]<t_max
        f_mu.=reaction_rate.*vec(prod(binomial.(state,sml_t);dims=1)) #计算反应倾向函数。比如某行mX+nY(+0*Z)->kM，X有a个，Y有b个，Z有c个，则计算C(a,m)*C(b,n)*C(c,0)*rate            
        F=sum(f_mu)
        r1=rand()
        tau=-log(r1)/F #产生反应时间
        r2=rand()
        f_mu_cumsum=cumsum(f_mu)
        next_mu=findfirst(x -> x>F*r2,f_mu_cumsum) #产生反应通道
        state.=state.-stoichiometry_matrix_left[next_mu,:].+stoichiometry_matrix_right[next_mu,:] #状态变化
        push!(state_record,state[:]) #记录状态
        push!(t,t[end]+tau) #记录时间节点
    end
    return (state_record,t)
end

#=
stoich_matrix_left=[ #消耗化学计量数
    #off on X
    1 0 0
    0 1 0
    0 1 0
    0 0 1
]
stoich_matrix_right=[ #生成化学计量数
    #off on X
    0 1 0
    1 0 0
    0 1 1
    0 0 0
]
rate=[0.5,0.3,0.8,0.1] #反应常数。注意用逗号间隔
t_max=1000000 #最大反应时间
state=[1,0,0] #初始各反应物状态。注意用逗号间隔

ss,tt=Gillespie(stoich_matrix_left,stoich_matrix_right,rate,state,t_max);
=#