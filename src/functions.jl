"""
    get_kernel(size=20, kernel="LPF")

Generate a filtering kernel
"""
function get_kernel(size=20, kernel="LPF")
    if size <= 0
        error("Invalid kernel size, it must be a positive integer")
    end
    if kernel == "LPF"
        kernel_vec = ones(size)  # flat averaging (no decay)
    elseif kernel == "COS"
        end_val = pi / 2.0
        kernel_vec = cos.(0.0 : end_val / size : end_val)  # cosine decay
    else
        error("Invalid kernel type, 'LPF' and 'COS' are implemented")
    end
    kernel_vec ./= sum(kernel_vec)  # normalize energy
    return kernel_vec
end

"""
    lic_process(field, kernel_size=20, kernel_type="LPF", seed=1, interpolate=true)

Returns a Line Integral Convolution image of the field
"""
function lic_process(field; kernel_size=20, kernel_type="LPF", seed=1, interpolate=true)

    image_width, image_height = size(field)

    # Prepare filtering kernel
    kernel = get_kernel(kernel_size, kernel_type)

    # Generate base random image
    random_image = rand(MersenneTwister(seed), image_width, image_height) .- 0.5
    seed += 1

    # Initialize LIC image
    img = zeros(image_width, image_height)

    # Compute LIC
    theta_mat = atan.(field[:, :, 1], field[:, :, 2])
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

    field_modulus = abs.(field[:,:,1] .^ 2 + field[:,:,2] .^ 2)
    field_log = 10.0*log10.(field_modulus)

    v, _ = findmax(field_log)
    vmin, _ = findmin(field_log)
    int_img = (field_log .- vmin) ./ (v - vmin)  # Intensity (normalize)

    alpha_lic = 0.3
    alpha_int = 1.0 - alpha_lic
    r = img .* alpha_lic .+ int_img .* alpha_int
    g = img .* alpha_lic .+ int_img .* alpha_int
    b = img .* alpha_lic .+ (1.0 .- int_img) .* alpha_int
    img1d = colorview(RGB, r, g, b)

    return img1d
end

"""
    lic(field, kernel_size=20, kernel_type="LPF", seed=1, interpolate=true)

Returns a Line Integral Convolution image of the field
"""
function lic(field; kernel_size=20, kernel_type="LPF", seed=1, interpolate=true)

    lic_img = lic_process(
        field, kernel_size=kernel_size, 
        kernel_type=kernel_type, seed=seed, 
        interpolate=interpolate
    )

    heatmap(lic_img)
end