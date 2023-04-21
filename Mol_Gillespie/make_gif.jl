#构建动画
export make_gif


function make_gif(frame_rec::Array{Float64, 3},path,fps,title="example",label="data")
    ani=Animation()
    x_min,x_max,y_min,y_max,z_min,z_max=frame_rec[1,1,1],frame_rec[1,1,1],frame_rec[1,2,1],frame_rec[1,2,1],frame_rec[1,3,1],frame_rec[1,3,1]
    for i in 1:size(frame_rec,3)
        x_min=min(x_min,minimum(frame_rec[:,1,i]))
        y_min=min(y_min,minimum(frame_rec[:,2,i]))
        z_min=min(z_min,minimum(frame_rec[:,3,i]))
        x_max=max(x_max,maximum(frame_rec[:,1,i]))
        y_max=max(y_max,maximum(frame_rec[:,2,i]))
        z_max=max(z_max,maximum(frame_rec[:,3,i]))
    end
    x_min=floor(Int64,x_min)
    y_min=floor(Int64,y_min)
    z_min=floor(Int64,z_min)
    x_max=ceil(Int64,x_max)
    y_max=ceil(Int64,y_max)
    z_max=ceil(Int64,z_max)
    for i in 1:size(frame_rec,3)
        p=plot(frame_rec[:,1,i],frame_rec[:,2,i],frame_rec[:,3,i],xlim=(x_min,x_max),ylim=(y_min,y_max),zlim=(z_min,z_max),title=title,label=false,lc=:black,lw=2,la=1);
        scatter!(frame_rec[:,1,i],frame_rec[:,2,i],frame_rec[:,3,i],label=label,legend=:outerright,mc=:red,ms=5,ma=1);
        frame(ani,p)
    end
    gif(ani,path,fps=fps)
end