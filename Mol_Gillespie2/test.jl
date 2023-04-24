#测试示例
using Plots
using LinearAlgebra
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
react_rate=[5,3,8,1]*1e-3 #反应常数。注意用逗号间隔
state=[1,0,0] #初始各反应物状态。注意用逗号间隔
distance_affect_rate=d->([max(0,5-2d),max(0,2-5d),0.,0.]*1e-3)

EP_index=[2,4]
part_num=10
mass=10000*ones(part_num)
pos=zeros(part_num,3);pos[1,:]=randn(3);pos[2,:]=pos[1,:]+normalize(randn(3))*(30+randn()); #第二个分子位置=前一个位置+随机方向*距离

offset_angle=3*pi/12 #用于描述三个单体形成的最大外角
for i in 2:part_num-1
    pos[i+1,:]=2*pos[i,:]-pos[i-1,:]+normalize(randn(3))*norm(pos[i,:]-pos[i-1,:])*sin(rand()*offset_angle) #r(i+1)=r(i)+r(i)-r(i-1)+normalize(randn(3))*norm(r(i)-r(i-1))sin(rand()*α)
end
initial_vel=randn(part_num,3)

M=zeros(part_num,part_num);M.=-1;for i in 1:part_num-1 M[i,i+1]=1e0 end
L0=30*ones(part_num,part_num)
#cq=sign(randn(part_num))
force_types=[Lennard_Jones(10,10),Rouse_model(M,L0)]
#force_types=[WCA(1e-1,20),Rouse_model(M,L0)]
#force_types=[Lennard_Jones(1e5,10)]
#force_types=[Lennard_Jones(1,10),FENE(M,30,11)]
#force_types=[Coulomb(cq),Lennard_Jones(1,1),Rouse_model(M,L0)]
#force_types=[Rouse_model(M,L0)]
t_max=convert(Int64,1e6)
boundary=0
iterate_method="langevin leap-frog"
#iterate_method="velocity verlet"

frame_interval=1000 #表示用于制作GIF的帧相邻间隔多少个单位时间，不需要制作GIF的话frame_interval就填Inf
#t_max/frame_interval最好不要大于1000，否则GIF制作会耗时半分钟以上



#能量最小化
#pos=min_energy(10000,1e-30,part_num,pos,mass,force_types,boundary);


#模拟
@time pos,state,frame_rec,state_record,t_record=molecular_simulation_and_Gillespie(state,stoich_matrix_left,stoich_matrix_right,react_rate,EP_index,distance_affect_rate,  t_max,part_num,mass,pos,initial_vel,force_types,iterate_method,boundary,frame_interval);

#制作gif
#make_gif(frame_rec,"./MDexample.gif",30,EP_index,"example","atom");