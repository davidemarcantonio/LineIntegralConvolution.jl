"""
    get_kernel(size=20, kernel="LPF")

Generate a filtering kernel
"""
function get_kernel(size=20, kernel="LPF")
    if kernel == "LPF"
        kernel_vec = ones(size)
    elseif kernel == "COS"
        kernel_vec = cos.(0.0 : pi / (2.0* size) : pi/2.0)
    else
        error("Invalid kernel")
    end
    kernel_vec ./= sum(kernel_vec)
    return kernel_vec
end

"""
    lic(field::Field2D, kernel_size=20, kernel_type="LPF", seed=1, interpolate=true)

Returns a Line Integral Convolution image of the field
"""
function lic_process(field::Field2D, kernel_size=20, kernel_type="LPF", seed=1, interpolate=true)

    image_width = length(field.pos_y)
    image_height = length(field.pos_x)
    f_lin = field.val

    # Prepare filtering kernel
    kernel = get_kernel(kernel_size, kernel_type)

    # Generate base random image
    random_image = rand(MersenneTwister(seed), image_width, image_height) .- 0.5
    seed += 1

    # Initialize LIC image
    img = zeros(image_width, image_height)

    # Compute LIC
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

"""
    lic(field::Field2D, kernel_size=20, kernel_type="LPF", seed=1, interpolate=true)

Returns a Line Integral Convolution image of the field
"""
function lic(field::Field2D, kernel_size=20, kernel_type="LPF", seed=1, interpolate=true)

    lic_img = lic_process(field, kernel_size, kernel_type, seed, interpolate)  # LIC

    # field_modulus = sum(abs.(values))
    # field_log = 10.0*log10.(field_modulus)

    v, i = findmax(field_result.z)
    vmin, i = findmin(field_result.z)
    int_img = (field_result.z .- vmin) ./ (v - vmin) # Intensity (normalize)

    alpha_lic = 0.3
    alpha_int = 1.0 - alpha_lic
    r = lic_img .* alpha_lic .+ int_img .* alpha_int
    g = lic_img .* alpha_lic .+ int_img .* alpha_int
    b = lic_img .* alpha_lic .+ (1.0 .- int_img) .* alpha_int
    img1d = colorview(RGB, r, g, b)

    heatmap(img1d)
end