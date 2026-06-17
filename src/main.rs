use std::{collections::HashMap, env, fs};
use sola_raylib::prelude::*;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        panic!("No file provided");
    }
    let file_path = &args[1];
    let contents = fs::read_to_string(file_path).expect("Unable to read file");

    let width = 1600;
    let height = 900;
    let global_scale = height as f64 / 16.0;
    let font_size = 20;

    let leech = generate_leech();

    let mut lines = contents.lines();
    let mut line = lines.next().expect("Empty file");

    // Obtain slicing vectors
    let mut slicing_vectors: Vec<[i32; 24]> = Vec::new();
    while !line.is_empty() {
        slicing_vectors.push(line.split_whitespace()
            .map(|x| x.parse::<i32>().expect("Slicing vector contains non-numbers"))
            .collect::<Vec<i32>>()
            .as_slice().try_into().expect("Invalid slicing vector length"));
        line = lines.next().expect("No projection vectors provided");
    }

    // Obtain projection vectors
    let projection_vectors: Vec<[i32; 24]> = lines.map(|line| line.split_whitespace()
            .map(|x| x.parse::<i32>().expect("Projection vector contains non-numbers"))
            .collect::<Vec<i32>>()
            .as_slice().try_into().expect("Invalid projection vector length"))
        .collect();
    if projection_vectors.len() < 2 {
        panic!("Not enough projection vectors");
    } else if projection_vectors.len() % 2 == 1 {
        panic!("Unpaired projection vector");
    }

    // Initial slicing
    let mut slice_heights = vec![0; slicing_vectors.len()];
    let mut active_slice_idx: usize = 0;
    let mut leech_slice: Vec<[i32; 24]> = slice(&leech, &slicing_vectors, &slice_heights);

    // Produce orthogonal basis vectors for slicing vector subspace
    let mut ortho_slicing_vectors: Vec<[f64; 24]> = slicing_vectors.iter().copied().map(|v| int_to_float_vector(v)).collect();
    if !ortho_slicing_vectors.is_empty() {
        let m = fdot(ortho_slicing_vectors[0], ortho_slicing_vectors[0]);
        for i in 0..24 {
            ortho_slicing_vectors[0][i] /= m;
        }
        for i in 0..ortho_slicing_vectors.len() {
            for j in i+1..ortho_slicing_vectors.len() {
                ortho_slicing_vectors[j] = vectors_to_basis(ortho_slicing_vectors[i], ortho_slicing_vectors[j]).0;
            }
        }
    }

    // Set projection vectors to be perpendicular to slicing vectors
    let disp_projection_vectors: Vec<[f64; 24]> = projection_vectors.iter().copied().map(|v| {
        let mut vf = int_to_float_vector(v);
        for s in ortho_slicing_vectors.iter().copied() {
            let v_proj = fdot(s, vf) / fdot(s, s);
            for i in 0..24 {
                vf[i] -= v_proj * s[i];
            }   
        }
        vf
    }).collect();
    let mut active_proj_idx = 0;

    // Count overlapping points in the projection and store representative
    let mut projected_points = HashMap::new();
    calculate_projected_points(&leech_slice, &projection_vectors, active_proj_idx, &mut projected_points);

    
    let (mut rl, thread) = sola_raylib::init()
        .size(width, height)
        .title("Leech Slicer")
        .build();

    let mut t:f32 = 1.0;
    let mut old_proj_idx = active_proj_idx;
    let switch_time = 0.5;

    while !rl.window_should_close() {

        t += rl.get_frame_time();

        // Projection switching
        if rl.is_key_pressed(KeyboardKey::KEY_D) || rl.is_key_pressed(KeyboardKey::KEY_A) {
            old_proj_idx = active_proj_idx;
            let new_idx = if rl.is_key_pressed(KeyboardKey::KEY_D) {active_proj_idx + 1} else {active_proj_idx + projection_vectors.len() - 1};
            active_proj_idx = new_idx.rem_euclid(projection_vectors.len() / 2);
            calculate_projected_points(&leech_slice, &projection_vectors, active_proj_idx, &mut projected_points);
            t = 0.0;
        }

        // Parallel slices
        if slice_heights.len() != 0 {

            // Change active slice vector
            if rl.is_key_pressed(KeyboardKey::KEY_DOWN) {
                active_slice_idx = (active_slice_idx + 1).rem_euclid(slice_heights.len());
            } else if rl.is_key_pressed(KeyboardKey::KEY_UP) {
                active_slice_idx = (active_slice_idx + slice_heights.len() - 1).rem_euclid(slice_heights.len());
            }

            // Increase slice height
            if rl.is_key_pressed(KeyboardKey::KEY_RIGHT) {
                let mut next_slice_height = i32::MAX;
                for v in &leech {
                    let v_height = dot(*v, slicing_vectors[active_slice_idx]);
                    if v_height > slice_heights[active_slice_idx] && v_height < next_slice_height {
                        next_slice_height = v_height;
                    }
                }
                if next_slice_height != i32::MAX {
                    slice_heights[active_slice_idx] = next_slice_height;
                    leech_slice = slice(&leech, &slicing_vectors, &slice_heights);
                    calculate_projected_points(&leech_slice, &projection_vectors, active_proj_idx, &mut projected_points);
                }
            }

            // Decrease slice height
            else if rl.is_key_pressed(KeyboardKey::KEY_LEFT) {
                let mut next_slice_height = i32::MIN;
                for v in &leech {
                    let v_height = dot(*v, slicing_vectors[active_slice_idx]);
                    if v_height < slice_heights[active_slice_idx] && v_height > next_slice_height {
                        next_slice_height = v_height;
                    }
                }
                if next_slice_height != i32::MIN {
                    slice_heights[active_slice_idx] = next_slice_height;
                    leech_slice = slice(&leech, &slicing_vectors, &slice_heights);
                    calculate_projected_points(&leech_slice, &projection_vectors, active_proj_idx, &mut projected_points);
                }
            }
        }

        let mut d = rl.begin_drawing(&thread);

        d.clear_background(Color::WHITE);
        
        if t as f64 >= switch_time {
            // Draw points with overlap counts
            let (b1, b2) = vectors_to_basis(disp_projection_vectors[2 * active_proj_idx], disp_projection_vectors[2 * active_proj_idx + 1]);
            for point in projected_points.values() {
                let count = point.0;
                let rep = int_to_float_vector(point.1);

                let x = fdot(b1, rep);
                let y = fdot(b2, rep);
                d.draw_circle(width/2 + (x * global_scale) as i32, height/2 - (y * global_scale) as i32, 5.0, Color::BLACK);

                d.draw_text(count.to_string().as_str(), width/2 + 5 + (x * global_scale) as i32, height/2 + 5 - (y * global_scale) as i32, font_size, Color::BLACK)
            }
        } else {
            // Projection switching animation
            // Draw all points separately, counting overlap is unreliable
            let t_clipped = f64::min(t as f64, switch_time);
            let mut disp_1 = [0.0; 24];
            let mut disp_2 = [0.0; 24];
            for i in 0..24 {
                disp_1[i] = t_clipped * disp_projection_vectors[2 * active_proj_idx][i] + (switch_time - t_clipped) * disp_projection_vectors[2 * old_proj_idx][i];
                disp_2[i] = t_clipped * disp_projection_vectors[2 * active_proj_idx + 1][i] + (switch_time - t_clipped) * disp_projection_vectors[2 * old_proj_idx + 1][i];
            }
            let (b1, b2) = vectors_to_basis(disp_1, disp_2);
            for point in leech_slice.iter().map(|v| int_to_float_vector(*v)) {
                let x = fdot(b1, point);
                let y = fdot(b2, point);
                d.draw_circle(width/2 + (x * global_scale) as i32, height/2 - (y * global_scale) as i32, 5.0, Color::BLACK);
            }
        }

        // Slice info
        for i in 0..slice_heights.len() {
            let height_str;
            if i == active_slice_idx {
                height_str = format!(">> {} <<", slice_heights[i]);
            } else {
                height_str = format!("   {}", slice_heights[i]);
            }
            d.draw_text(height_str.as_str(), font_size, 2 * font_size * (i as i32 + 1), font_size, Color::BLACK);
        }
    }
}

fn generate_leech() -> Vec<[i32; 24]> {
    let mut leech = Vec::with_capacity(196560);

    let golay = generate_golay();

    // Vectors of type (4^2, 0^22)
    for i in 0..24 {
        for j in i+1..24 {
            let mut v = [0; 24];
            v[i] = 4;
            v[j] = 4;
            leech.push(v);
            v[j] = -4;
            leech.push(v);
            v[i] = -4;
            leech.push(v);
            v[j] = 4;
            leech.push(v);
        }
    }

    // Vectors of type (2^8, 0^16)
    for g in &golay {
        if g.iter().sum::<i32>() == 8 {
            for i in 0..256 {
                let mut v = *g;
                let mut counter = 7;
                for idx in 0..24 {
                    if v[idx] != 0 {
                        if i & (1 << counter) == 0 {
                            v[idx] *= 2;
                        } else {
                            v[idx] *= -2;
                        }
                        counter -= 1;
                    }
                }
                if v.iter().sum::<i32>().rem_euclid(8) == 0 {
                    leech.push(v);
                }
            }
        }
    }

    // Vectors of type (-3, 1^23)
    for i in 0..24 {
        let mut v_init = [1; 24];
        v_init[i] = -3;
        for g in &golay {
            let mut v = v_init;
            for idx in 0..24 {
                if g[idx] == 1 {
                    v[idx] *= -1;
                }
            }
            leech.push(v);
        }
    }

    leech
}

fn generate_golay() -> Vec<[i32; 24]> {
    let golay_matrix = [
        [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],
        [1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0],
        [1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0],
        [0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1],
    ];

    (0..4096).map(|x| {
        let mut c = [0; 24];
        for n in 0..12 {
            if x & (1 << n) != 0 {
                for i in 0..24 {
                    c[i] ^= golay_matrix[n][i];
                }
            }
        }
        c
    }).collect()
}

fn slice(points: &Vec<[i32; 24]>, slicing_vectors: &Vec<[i32; 24]>, slice_heights: &Vec<i32>) -> Vec<[i32; 24]> {
    points.iter().copied().filter(|x| {
        for i in 0..slicing_vectors.len() {
            if dot(*x, slicing_vectors[i]) != slice_heights[i] {
                return false;
            }
        }
        true
    }).collect()
}

fn calculate_projected_points(points: &Vec<[i32; 24]>, projection_vectors: &Vec<[i32; 24]>, projection_index: usize, projected_points: &mut HashMap<(i32, i32), (i32, [i32; 24])>) {
    projected_points.clear();
    for i in 0..points.len() {
        projected_points.entry((dot(points[i], projection_vectors[2 * projection_index]), dot(points[i], projection_vectors[2 * projection_index + 1])))
            .and_modify(|c: &mut (i32, [i32; 24])| {c.0 += 1})
            .or_insert((1, points[i]));
    }
}

fn dot(a: [i32; 24], b: [i32; 24]) -> i32 {
    let mut result = 0;
    for i in 0..24 {
        result += a[i] * b[i];
    }
    result
}

fn fdot(a: [f64; 24], b: [f64; 24]) -> f64 {
    let mut result = 0.0;
    for i in 0..24 {
        result += a[i] * b[i];
    }
    result
}

// okay yes this is a bad name for this function
// I'm sorry </3
fn vectors_to_basis(v1: [f64; 24], v2: [f64; 24]) -> ([f64; 24], [f64; 24]) {
    let v2_proj = fdot(v1, v2) / fdot(v1, v1);
    let mut vx = [0.0; 24];
    let mut vy = v1;
    for i in 0..24 {
        vx[i] = v2[i] - vy[i] * v2_proj;
    }
    let mx = f64::sqrt(fdot(vx, vx));
    let my = f64::sqrt(fdot(vy, vy));
    for i in 0..24 {
        vx[i] /= mx;
        vy[i] /= my;
    }
    (vx, vy)
}

fn int_to_float_vector(v: [i32; 24]) -> [f64; 24] {
    let mut vf = [0.0; 24];
    for i in 0..24 {
        vf[i] = v[i] as f64;
    }
    vf
}