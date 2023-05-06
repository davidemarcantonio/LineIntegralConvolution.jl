using Random

"""
    greet_your_package_name()

test
"""
function greet_your_package_name()
    return "Hello world!"
end

"""
    meshgrid(x_range, y_range)

Returns a 2D grid of points.
"""
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

"""
    simulate_2d_electrostatic(
        x_min, x_max, y_min, y_max, pt_per_meter, n_charges=1, 
        distribution="circle", seed=1, charge_value=1e-6)

Returns a 2D grid of points.
"""
function simulate_2d_electrostatic(
    x_min, x_max, y_min, y_max, pt_per_meter, n_charges=1, 
    distribution="circle", seed=1, charge_value=1e-6)

    epsilon_zero = 8.85e-12  # empty space dielectric constant
    k = 1.0 / (4.0 * pi * epsilon_zero)  # constant, see electric field generated form a point source

    step_m = 1.0 / pt_per_meter
    pts_x = x_min : step_m : x_max
    pts_y = y_min : step_m : y_max

    # num_x_points = (x_max - x_min) * pt_per_meter + 1  # +1 to include both ends
    # num_y_points = (y_max - y_min) * pt_per_meter + 1  # +1 to include both ends
    num_x_points = length(pts_x)
    num_y_points = length(pts_y)

    points = meshgrid(pts_x, pts_y)
    # points.x
    # points.y

    field_lin = zeros(num_y_points, num_x_points, 2)  # initialize zero electric field 
    field_modulus = zeros(num_y_points, num_x_points)  # initialize the modulus only image

    if distribution == "random"
        # generate a random position [x, y] for the new charges
        # the code below forces the position to be 'in the center'
        charge_x = rand(MersenneTwister(seed), Float64, n_charges) * (x_max-x_min) .+ x_min
        seed += 1
        charge_y = rand(MersenneTwister(seed), Float64, n_charges) * (y_max-y_min) .+ y_min
        seed += 1
        charge_q = rand(MersenneTwister(seed), (-1, 1), n_charges) .* charge_value
        seed += 1
    elseif distribution == "circle"
        if n_charges > 1
            angles = 0.0 : 2.0pi / n_charges : 1.999pi
            c_cos = cos.(angles)
            c_sin = sin.(angles)
        else
            c_cos = 0.0
            c_sin = 0.0
        end
        charge_x = x_min + (x_max-x_min) / 2 .+ c_cos * (x_max-x_min) / 4 .+ rand(MersenneTwister(seed), Float64) * step_m / 100
        charge_y = y_min + (y_max-y_min) / 2 .+ c_sin * (y_max-y_min) / 4 .+ rand(MersenneTwister(seed), Float64) * step_m / 100
        charge_q = ones(n_charges) * charge_value
    else
        error("Not implemented!")
    end

    for i = 1 : n_charges  # for each electric charge
        # for each pixel of the image, compute the electric field in that point
        # relatively to the currently generated charge
    
        # compute the distance of the current pixel position with
        # respect to the charge position, scaling to make them 'physical'
        dist_x = points.x .- charge_x[i]
        dist_y = points.y .- charge_y[i]
        # compute the squared linear distance between source and pixel
        dist2 = dist_x .^ 2 .+ dist_y .^ 2

        # compute the electric field value for the current pixel
        values = (k * charge_value) ./ dist2
        
        field_modulus .+= abs.(values)

        den = abs.(dist_x) .+ abs.(dist_y)
        x_versor = dist_x ./ den
        y_versor = dist_y ./ den

        # add the field computed to the overall field in this point 
        # (superposition principle)
        field_lin[:, :, 1] .+= values .* x_versor
        field_lin[:, :, 2] .+= values .* y_versor
    end

    field_log = 10.0*log10.(field_modulus)

    return (x=pts_x, y=pts_y, z=field_log, f=field_lin, cx=charge_x, cy=charge_y, cq=charge_q)
end

"""
Field2D

"""
struct Field2D
    pos_x
    pos_y
    val_x
    val_y
end

"""
    lic(f_x, f_y, f_lin, kernel_size=20, kernel_type="LPF", seed=1, interpolate=true)

lic(f_x, f_y, f_lin, kernel_size=20, kernel_type="LPF", seed=1, interpolate=true)

Returns a 2D grid of points.
"""
function lic(f_x, f_y, f_lin, kernel_size=20, kernel_type="LPF", seed=1, interpolate=true)

    if kernel_type == "LPF"
        kernel = ones(kernel_size)
    elseif kernel_type == "COS"
        kernel = cos.(0.0 : pi / (2.0* kernel_size) : pi/2.0)
    else
        error("Invalid kernel")
    end

    kernel ./= sum(kernel)

    image_width = length(f_y)
    image_height = length(f_x)

    ## initial random image
    random_image = rand(MersenneTwister(seed), image_width, image_height) .- 0.5
    seed += 1

    # new images
    img = zeros(image_width, image_height)

    # compute LIC

    theta_mat = atan.(f_lin[:, :, 1], f_lin[:, :, 2])
    Δx = sin.(theta_mat)
    Δy = cos.(theta_mat)
    if interpolate
        abs_sin = abs.(Δx)
        abs_cos = abs.(Δy)
    end

    target_num_prints = 10

    for x = 1 : image_width
        for y = 1 : image_height

            if rand() > 1.0 - target_num_prints / (image_width * image_height)
                print((((x * (image_height-1) + y)/(image_width * image_height))*100), " %\n");
            end

            c_pix = random_image[x, y] * kernel[1]
            xp = x
            yp = y
            cdx = 0
            cdy = 0

            for w = 2 : kernel_size

                cdx += Δx[xp, yp]
                cdy += Δy[xp, yp]

                # compute new position pixel
                xpc = x + Int(round(cdx))
                ypc = y + Int(round(cdy))

                # check new position is valid
                if xpc < 1
                    xpc = 1
                elseif xpc > image_width
                    xpc = image_width
                end
                if ypc < 1
                    ypc = 1
                elseif ypc > image_height
                    ypc = image_height
                end

                if interpolate
                    norm_factor = (1 - abs_cos[xp, yp]) * (1 - abs_sin[xp, yp])
                    c_pix += random_image[xp, yp] * kernel[w] * norm_factor
                    norm_factor = abs_cos[xp, yp]
                    c_pix += random_image[xpc, yp] * kernel[w] * norm_factor
                    norm_factor = abs_sin[xp, yp]
                    c_pix += random_image[xp, ypc] * kernel[w] * norm_factor
                    norm_factor = abs_cos[xp, yp] * abs_sin[xp, yp]
                    c_pix += random_image[xpc, ypc] * kernel[w] * norm_factor
                else
                    c_pix += random_image[xpc, ypc] * kernel[w]
                end

                xp = xpc
                yp = ypc
            end

            img[x, y] = c_pix
        end
    end

    img .+= 0.5

    return (rnd_img=random_image, final_img=img)
end