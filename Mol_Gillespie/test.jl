#测试示例

using Random
Random.seed!(123)
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
react_rate=[0.5,0.3,0.8,0.1] #反应常数。注意用逗号间隔
initial_state=[1,0,0] #初始各反应物状态。注意用逗号间隔
distance_affect_rate=d->[max(0.,3-d),max(0.,10-d),max(0.,6-2d),max(0.,7-0.5d)]

EP_index=[2,4]
part_num=10
mass=ones(part_num)
initial_pos=zeros(part_num,3);initial_pos[1,:]=2*randn(3);
for i in 2:part_num
    initial_pos[i,:]=initial_pos[i-1,:]+2*randn(3)
end
initial_vel=1*randn(part_num,3)

M=zeros(part_num,part_num);M.=-1;for i in 1:part_num-1 M[i,i+1]=5*rand() end
L0=ones(part_num,part_num)
force_types=[Lennard_Jones(0.25,0.5),Rouse_model(M,L0)]
#force_types=[Lennard_Jones(0.25,0.1),FENE(M,10,1)]
dt=1e-3
t_max=150
boundary=0
iterate_method="langevin leap-frog"
#iterate_method="velocity verlet"

frame_interval=100 #表示用于制作GIF的帧相邻间隔多少个dt，不需要制作GIF的话frame_interval就填Inf
#t_max/(dt*frame_interval)最好不要大于1000，否则GIF制作会耗时半分钟以上



#能量最小化
pos=min_energy(10000,100,part_num,initial_pos,mass,force_types,boundary);

#模拟
@time pos,state,frame_rec,state_record,t_record=molecular_simulation_and_Gillespie(initial_state,stoich_matrix_left,stoich_matrix_right,react_rate,EP_index,distance_affect_rate,  t_max,dt,part_num,mass,pos,initial_vel,force_types,iterate_method,boundary,frame_interval);
#制作gif
#make_gif(frame_rec,"./MDexample.gif",30,"example","atom");