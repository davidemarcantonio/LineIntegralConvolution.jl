using Random

function greet_your_package_name()
    return "Hello world!"
end

function meshgrid(x_range, y_range)
    num_x_pts = length(x_range)
    num_y_pts = length(y_range)
    pts_x = zeros(num_y_pts, num_x_pts)
    pts_y = zeros(num_y_pts, num_x_pts)
    for cnt_x = 1 : num_x_pts
        for cnt_y = 1 : num_y_pts
            pts_x[cnt_y, cnt_x] = x_range[cnt_x]
            pts_y[cnt_y, cnt_x] = y_range[cnt_y]
        end
    end
    return (x=pts_x, y=pts_y)
end