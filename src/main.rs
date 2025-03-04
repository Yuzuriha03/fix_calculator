use std::f64::consts::PI;
use eframe::egui;

const A: f64 = 6378137.0;
const B: f64 = 6356752.314245;
const F: f64 = 1.0 / 298.257223563;
const TOLERANCE: f64 = 1e-12;

// Simple conversion between degrees and radians
fn deg_to_rad(deg: f64) -> f64 {
    deg * PI / 180.0
}
fn rad_to_deg(rad: f64) -> f64 {
    rad * 180.0 / PI
}

/// Vincenty's Inverse Method: Given two points, compute distance (m) and initial azimuth (°).
fn vincenty_inverse(lat1: f64, lon1: f64, lat2: f64, lon2: f64) -> Result<GeodesicInverseResult, String> {
    // Convert degrees to radians
    let phi1 = deg_to_rad(lat1);
    let phi2 = deg_to_rad(lat2);
    let lambda1 = deg_to_rad(lon1);
    let lambda2 = deg_to_rad(lon2);

    let u1 = ((1.0 - F) * phi1.tan()).atan();
    let u2 = ((1.0 - F) * phi2.tan()).atan();
    let sin_u1 = u1.sin();
    let cos_u1 = u1.cos();
    let sin_u2 = u2.sin();
    let cos_u2 = u2.cos();

    let l = lambda2 - lambda1;
    let mut lambda = l;
    let mut iter_count = 0;
    let (mut sin_sigma, mut cos_sigma, mut sigma, mut sin_alpha, mut cos_sq_alpha, mut cos2sigma_m, mut c);
    loop {
        iter_count += 1;
        if iter_count > 1000 {
            return Err("Vincenty's inverse formula did not converge".to_string());
        }
        let sin_lambda = lambda.sin();
        let cos_lambda = lambda.cos();
        let temp1 = cos_u2 * sin_lambda;
        let temp2 = cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda;
        sin_sigma = (temp1 * temp1 + temp2 * temp2).sqrt();
        if sin_sigma == 0.0 {
            return Ok(GeodesicInverseResult { s12: 0.0, azi1: 0.0 });
        }
        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lambda;
        sigma = sin_sigma.atan2(cos_sigma);
        sin_alpha = cos_u1 * cos_u2 * sin_lambda / sin_sigma;
        cos_sq_alpha = 1.0 - sin_alpha * sin_alpha;
        cos2sigma_m = if cos_sq_alpha != 0.0 {
            cos_sigma - 2.0 * sin_u1 * sin_u2 / cos_sq_alpha
        } else {
            0.0
        };
        c = F / 16.0 * cos_sq_alpha * (4.0 + F * (4.0 - 3.0 * cos_sq_alpha));
        let lambda_prev = lambda;
        lambda = l + (1.0 - c) * F * sin_alpha * (sigma + c * sin_sigma * (cos2sigma_m + c * cos_sigma * (-1.0 + 2.0 * cos2sigma_m * cos2sigma_m)));
        if (lambda - lambda_prev).abs() < TOLERANCE {
            break;
        }
    }
    let u_sq = cos_sq_alpha * (A * A - B * B) / (B * B);
    let a_coef = 1.0 + u_sq / 16384.0 * (4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0 * u_sq)));
    let b_coef = u_sq / 1024.0 * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)));
    let delta_sigma = b_coef * sin_sigma * (cos2sigma_m + b_coef / 4.0 * (cos_sigma * (-1.0 + 2.0 * cos2sigma_m * cos2sigma_m)
        - b_coef / 6.0 * cos2sigma_m * (-3.0 + 4.0 * sin_sigma * sin_sigma) * (-3.0 + 4.0 * cos2sigma_m * cos2sigma_m)));
    let s = B * a_coef * (sigma - delta_sigma);
    let azi1 = (cos_u2 * lambda.sin()).atan2(cos_u1 * sin_u2 - sin_u1 * cos_u2 * lambda.cos());
    Ok(GeodesicInverseResult { s12: s, azi1: rad_to_deg(azi1) })
}

struct GeodesicInverseResult {
    s12: f64,  // Distance (m)
    azi1: f64, // Initial azimuth (°)
}

/// Vincenty's Direct Method: Given a start point, an initial azimuth (°) and distance (m), compute endpoint coordinates (°).
fn vincenty_direct(lat1: f64, lon1: f64, azi1: f64, s: f64) -> Result<(f64, f64), String> {
    let phi1 = deg_to_rad(lat1);
    let lambda1 = deg_to_rad(lon1);
    let alpha1 = deg_to_rad(azi1);

    let u1 = ((1.0 - F) * phi1.tan()).atan();
    let sin_u1 = u1.sin();
    let cos_u1 = u1.cos();
    let sigma1 = (u1.tan()).atan2(alpha1.cos());
    let sin_alpha = cos_u1 * alpha1.sin();
    let cos_sq_alpha = 1.0 - sin_alpha * sin_alpha;
    let u_sq = cos_sq_alpha * (A * A - B * B) / (B * B);
    let a_coef = 1.0 + u_sq / 16384.0 * (4096.0 + u_sq * (-768.0 + u_sq * (320.0 - 175.0 * u_sq)));
    let b_coef = u_sq / 1024.0 * (256.0 + u_sq * (-128.0 + u_sq * (74.0 - 47.0 * u_sq)));
    let mut sigma = s / (B * a_coef);
    let mut sigma_prev = 0.0;
    let mut iter_count = 0;
    let mut cos2sigma_m = 0.0;
    let mut sin_sigma = 0.0;
    let mut cos_sigma = 0.0;
    while (sigma - sigma_prev).abs() > TOLERANCE {
        iter_count += 1;
        if iter_count > 1000 {
            return Err("Vincenty's direct formula did not converge".to_string());
        }
        sigma_prev = sigma;
        sin_sigma = sigma.sin();
        cos_sigma = sigma.cos();
        cos2sigma_m = (2.0 * sigma1 + sigma).cos();
        let delta_sigma = b_coef * sin_sigma * (cos2sigma_m + b_coef / 4.0 * (cos_sigma * (-1.0 + 2.0 * cos2sigma_m * cos2sigma_m)
            - b_coef / 6.0 * cos2sigma_m * (-3.0 + 4.0 * sin_sigma * sin_sigma) * (-3.0 + 4.0 * cos2sigma_m * cos2sigma_m)));
        sigma = s / (B * a_coef) + delta_sigma;
    }
    let tmp = sin_u1 * cos_sigma + cos_u1 * sin_sigma * alpha1.cos();
    let phi2 = (tmp).atan2((1.0 - F) * ((sin_alpha * sin_alpha)
        + (sin_u1 * sin_sigma - cos_u1 * cos_sigma * alpha1.cos()).powi(2)).sqrt());
    let lambda = (sin_sigma * alpha1.sin()).atan2(cos_u1 * cos_sigma - sin_u1 * sin_sigma * alpha1.cos());
    let c = F / 16.0 * cos_sq_alpha * (4.0 + F * (4.0 - 3.0 * cos_sq_alpha));
    let l = lambda - (1.0 - c) * F * sin_alpha * (sigma + c * sin_sigma * (cos2sigma_m + c * cos_sigma * (-1.0 + 2.0 * cos2sigma_m * cos2sigma_m)));
    let lon2 = lambda1 + l;
    Ok((rad_to_deg(phi2), rad_to_deg(lon2)))
}

/// Golden Section Search: Finds the t in [a, b] that minimizes function f.
fn golden_section_search<F: Fn(f64) -> f64>(f: F, mut a: f64, mut b: f64, tol: f64) -> f64 {
    let gr = (5.0_f64.sqrt() + 1.0) / 2.0;
    let mut c = b - (b - a) / gr;
    let mut d = a + (b - a) / gr;
    while (b - a) > tol {
        if f(c) < f(d) {
            b = d;
        } else {
            a = c;
        }
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
    }
    (a + b) / 2.0
}

/// Parse a coordinate string supporting directional letters, DMS and compact formats.
fn parse_coordinate(input: &str, is_longitude: bool) -> Result<f64, String> {
    let input = input.trim().to_uppercase();
    if input.is_empty() {
        return Err("Coordinate input is empty".to_string());
    }
    // 如果以字母开头，则认为首字符为方向
    let first = input.chars().next().unwrap();
    if "NSEW".contains(first) {
        let direction = first;
        let remaining = input[1..].trim();
        let parts: Vec<&str> = remaining.split(|c| "°'\"’” ".contains(c))
            .filter(|s| !s.is_empty())
            .collect();
        if parts.len() > 1 {
            if let Ok(deg) = parts[0].parse::<f64>() {
                let minutes = if parts.len() >= 2 { parts[1].parse::<f64>().unwrap_or(0.0) } else { 0.0 };
                let seconds = if parts.len() >= 3 { parts[2].parse::<f64>().unwrap_or(0.0) } else { 0.0 };
                if !(0.0..60.0).contains(&minutes) || !(0.0..60.0).contains(&seconds) {
                    return Err("Minutes or seconds value out of range (0-60)".to_string());
                }
                let mut decimal = deg + minutes / 60.0 + seconds / 3600.0;
                if direction == 'S' || direction == 'W' {
                    decimal = -decimal;
                }
                if is_longitude && (decimal < -180.0 || decimal > 180.0) {
                    return Err("Longitude is out of valid range (-180° ~ 180°)".to_string());
                } else if !is_longitude && (decimal < -90.0 || decimal > 90.0) {
                    return Err("Latitude is out of valid range (-90° ~ 90°)".to_string());
                }
                return Ok(decimal);
            }
        }
        // 紧凑格式处理，例如 N394022.1
        let num_str: String = remaining.chars().filter(|c| c.is_digit(10) || *c == '.').collect();
        if num_str.is_empty() {
            return Err("The numeric part of the coordinate is empty".to_string());
        }
        let degree_digits = if is_longitude { 3 } else { 2 };
        if num_str.len() < degree_digits + 2 {
            return Err(format!("The numeric part of the coordinate is too short, requires at least {} digits for degrees plus 2 digits for minutes", degree_digits));
        }
        let deg: f64 = num_str[..degree_digits].parse().map_err(|_| "Invalid compact format numeric value".to_string())?;
        let minutes: f64 = num_str[degree_digits..degree_digits+2].parse().map_err(|_| "Invalid compact format numeric value".to_string())?;
        let seconds: f64 = if num_str.len() > degree_digits+2 {
            num_str[degree_digits+2..].parse().unwrap_or(0.0)
        } else { 0.0 };
        if !(0.0..60.0).contains(&minutes) || !(0.0..60.0).contains(&seconds) {
            return Err("Minutes or seconds value out of range (0-60)".to_string());
        }
        let mut decimal = deg + minutes / 60.0 + seconds / 3600.0;
        if direction == 'S' || direction == 'W' {
            decimal = -decimal;
        }
        if is_longitude && (decimal < -180.0 || decimal > 180.0) {
            return Err("Longitude is out of valid range (-180° ~ 180°)".to_string());
        } else if !is_longitude && (decimal < -90.0 || decimal > 90.0) {
            return Err("Latitude is out of valid range (-90° ~ 90°)".to_string());
        }
        return Ok(decimal);
    }
    // 尝试直接解析为浮点数
    if let Ok(val) = input.parse::<f64>() {
        return Ok(val);
    }
    // 拆分格式，例如 "39°40'22.1\"S" 或 "-39 40 22.1"
    let mut parts: Vec<&str> = input.split(|c| "°'\"’” ".contains(c))
        .filter(|s| !s.is_empty())
        .collect();
    let mut sign = 1.0;
    if parts.len() > 0 && (parts[0] == "-" || parts[0] == "+") {
        if parts[0] == "-" {
            sign = -1.0;
        }
        parts.remove(0);
    }
    let mut direction: Option<char> = None;
    if let Some(last) = parts.last() {
        if last.len() == 1 && "NSEW".contains(*last) {
            direction = Some(last.chars().next().unwrap());
            parts.pop();
        }
    }
    if let Some(dir) = direction {
        if dir == 'S' || dir == 'W' {
            sign = -1.0;
        }
    }
    let deg: f64 = parts.get(0).ok_or("Invalid coordinate component".to_string())?
        .parse().map_err(|_| "Invalid coordinate component".to_string())?;
    let minutes: f64 = if parts.len() >= 2 { parts[1].parse().unwrap_or(0.0) } else { 0.0 };
    let seconds: f64 = if parts.len() >= 3 { parts[2].parse().unwrap_or(0.0) } else { 0.0 };
    if !(0.0..60.0).contains(&minutes) || !(0.0..60.0).contains(&seconds) {
        return Err("Minutes or seconds value out of range (0-60)".to_string());
    }
    let decimal = sign * (deg + minutes / 60.0 + seconds / 3600.0);
    if is_longitude && (decimal < -180.0 || decimal > 180.0) {
        return Err("Longitude is out of valid range (-180° ~ 180°)".to_string());
    } else if !is_longitude && (decimal < -90.0 || decimal > 90.0) {
        return Err("Latitude is out of valid range (-90° ~ 90°)".to_string());
    }
    Ok(decimal)
}

/// Convert a decimal coordinate to a DMS string with directional letter.
fn format_dms(decimal: f64, is_lon: bool) -> String {
    let abs_val = decimal.abs();
    let deg = abs_val.floor() as i64;
    let rem = (abs_val - deg as f64) * 60.0;
    let minutes = rem.floor() as i64;
    let seconds = (rem - minutes as f64) * 60.0;
    let dir = if is_lon {
        if decimal >= 0.0 { "E" } else { "W" }
    } else {
        if decimal >= 0.0 { "N" } else { "S" }
    };
    format!("{}°{:02}'{:05.2}\"{}", deg, minutes, seconds, dir)
}

/// Core calculation logic: Compute the FI fix point based on runway and DME inputs.
fn calculate_fix(
    op_runway: (f64, f64),
    sm_runway: (f64, f64),
    dme1: (f64, f64),
    dme2: (f64, f64),
    dist1_m: f64,
    dist2_m: f64,
) -> Result<(f64, f64), String> {
    let inv = vincenty_inverse(op_runway.0, op_runway.1, sm_runway.0, sm_runway.1)?;
    let azimuth = inv.azi1;
    let base_dist = inv.s12;
    let get_point = |t: f64| -> Result<(f64, f64), String> {
        vincenty_direct(op_runway.0, op_runway.1, azimuth, t * base_dist)
    };
    let error = |t: f64| -> f64 {
        if let Ok((lat, lon)) = get_point(t) {
            let d1 = vincenty_inverse(lat, lon, dme1.0, dme1.1).map(|res| res.s12).unwrap_or(f64::INFINITY);
            let d2 = vincenty_inverse(lat, lon, dme2.0, dme2.1).map(|res| res.s12).unwrap_or(f64::INFINITY);
            (d1 - dist1_m).powi(2) + (d2 - dist2_m).powi(2)
        } else {
            f64::INFINITY
        }
    };
    let optimal_t = golden_section_search(error, 1.0, 10.0, 1e-6);
    get_point(optimal_t)
}

/// 以下辅助函数用于过滤输入内容
/// 对于Distance输入框：只允许数字和小数点
fn filter_distance(input: &str) -> String {
    input.chars().filter(|c| c.is_digit(10) || *c == '.').collect()
}
/// 对于经度输入框：字母部分只允许 E 和 W，其它字母过滤，其他字符保留
fn filter_longitude(input: &str) -> String {
    input.chars().filter(|c| {
        if c.is_alphabetic() {
            let uc = c.to_ascii_uppercase();
            uc == 'E' || uc == 'W'
        } else {
            true
        }
    }).collect()
}
/// 对于纬度输入框：字母部分只允许 N 和 S，其它字母过滤，其他字符保留
fn filter_latitude(input: &str) -> String {
    input.chars().filter(|c| {
        if c.is_alphabetic() {
            let uc = c.to_ascii_uppercase();
            uc == 'N' || uc == 'S'
        } else {
            true
        }
    }).collect()
}

/// Main application implementing the GUI with eframe/egui
struct AviationCoordinateCalculatorApp {
    // Input Parameter
    opposite_runway_lat: String,
    opposite_runway_lon: String,
    runway_lat: String,
    runway_lon: String,
    dme1_lat: String,
    dme1_lon: String,
    dme2_lat: String,
    dme2_lon: String,
    dist_dme1: String, // NM
    dist_dme2: String, // NM
    result: String,
    error: Option<String>,
}

impl Default for AviationCoordinateCalculatorApp {
    fn default() -> Self {
        Self {
            opposite_runway_lat: "".to_owned(),
            opposite_runway_lon: "".to_owned(),
            runway_lat: "".to_owned(),
            runway_lon: "".to_owned(),
            dme1_lat: "".to_owned(),
            dme1_lon: "".to_owned(),
            dme2_lat: "".to_owned(),
            dme2_lon: "".to_owned(),
            dist_dme1: "0.0".to_owned(),
            dist_dme2: "0.0".to_owned(),
            result: "".to_owned(),
            error: None,
        }
    }
}

impl eframe::App for AviationCoordinateCalculatorApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Input Filtering
        self.dist_dme1 = filter_distance(&self.dist_dme1);
        self.dist_dme2 = filter_distance(&self.dist_dme2);
        self.opposite_runway_lat = filter_latitude(&self.opposite_runway_lat);
        self.runway_lat = filter_latitude(&self.runway_lat);
        self.dme1_lat = filter_latitude(&self.dme1_lat);
        self.dme2_lat = filter_latitude(&self.dme2_lat);
        self.opposite_runway_lon = filter_longitude(&self.opposite_runway_lon);
        self.runway_lon = filter_longitude(&self.runway_lon);
        self.dme1_lon = filter_longitude(&self.dme1_lon);
        self.dme2_lon = filter_longitude(&self.dme2_lon);

        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.heading("Parameter Configuration");

                ui.separator();

                ui.group(|ui| {
                    ui.heading("Runway Coordinates");
                    egui::Grid::new("runway_grid")
                        .num_columns(2)
                        .spacing([40.0, 8.0])
                        .striped(true)
                        .show(ui, |ui| {
                            ui.label("Opposite Runway Latitude:");
                            ui.text_edit_singleline(&mut self.opposite_runway_lat);
                            ui.end_row();
                            ui.label("Opposite Runway Longitude:");
                            ui.text_edit_singleline(&mut self.opposite_runway_lon);
                            ui.end_row();
                            ui.label("Runway Latitude:");
                            ui.text_edit_singleline(&mut self.runway_lat);
                            ui.end_row();
                            ui.label("Runway Longitude:");
                            ui.text_edit_singleline(&mut self.runway_lon);
                            ui.end_row();
                        });
                });

                ui.separator();

                ui.group(|ui| {
                    ui.heading("DME Coordinates");
                    egui::Grid::new("dme_grid")
                        .num_columns(2)
                        .spacing([40.0, 8.0])
                        .striped(true)
                        .show(ui, |ui| {
                            ui.label("DME 1 Latitude:");
                            ui.text_edit_singleline(&mut self.dme1_lat);
                            ui.end_row();
                            ui.label("DME 1 Longitude:");
                            ui.text_edit_singleline(&mut self.dme1_lon);
                            ui.end_row();
                            ui.label("DME 2 Latitude:");
                            ui.text_edit_singleline(&mut self.dme2_lat);
                            ui.end_row();
                            ui.label("DME 2 Longitude:");
                            ui.text_edit_singleline(&mut self.dme2_lon);
                            ui.end_row();
                        });
                });

                ui.separator();

                ui.group(|ui| {
                    ui.heading("DME Distances (NM)");
                    egui::Grid::new("distance_grid")
                        .num_columns(2)
                        .spacing([40.0, 8.0])
                        .striped(true)
                        .show(ui, |ui| {
                            ui.label("Distance to DME 1:");
                            ui.text_edit_singleline(&mut self.dist_dme1);
                            ui.end_row();
                            ui.label("Distance to DME 2:");
                            ui.text_edit_singleline(&mut self.dist_dme2);
                            ui.end_row();
                        });
                });

                ui.separator();

                if ui.button("Calculate").clicked() {
                    self.error = None;
                    self.result.clear();

                    let parse = |s: &str, is_lon: bool| parse_coordinate(s, is_lon);
                    let op_runway = (
                        parse(&self.opposite_runway_lat, false),
                        parse(&self.opposite_runway_lon, true)
                    );
                    let sm_runway = (
                        parse(&self.runway_lat, false),
                        parse(&self.runway_lon, true)
                    );
                    let dme1 = (
                        parse(&self.dme1_lat, false),
                        parse(&self.dme1_lon, true)
                    );
                    let dme2 = (
                        parse(&self.dme2_lat, false),
                        parse(&self.dme2_lon, true)
                    );
                    let dist_dme1: Result<f64, _> = self.dist_dme1.trim().parse();
                    let dist_dme2: Result<f64, _> = self.dist_dme2.trim().parse();

                    if let (Ok(op_lat), Ok(op_lon), Ok(sm_lat), Ok(sm_lon),
                            Ok(d1_lat), Ok(d1_lon), Ok(d2_lat), Ok(d2_lon),
                            Ok(dd1), Ok(dd2)) =
                        (op_runway.0, op_runway.1, sm_runway.0, sm_runway.1,
                         dme1.0, dme1.1, dme2.0, dme2.1, dist_dme1, dist_dme2)
                    {
                        // Convert NM to meters (1 NM = 1852 m)
                        let dd1_m = dd1 * 1852.0;
                        let dd2_m = dd2 * 1852.0;
                        match calculate_fix((op_lat, op_lon), (sm_lat, sm_lon),
                                              (d1_lat, d1_lon), (d2_lat, d2_lon),
                                              dd1_m, dd2_m) {
                            Ok((fix_lat, fix_lon)) => {
                                let dms_lat = format_dms(fix_lat, false);
                                let dms_lon = format_dms(fix_lon, true);
                                let dec_lat = format!("{:.8}", fix_lat);
                                let dec_lon = format!("{:.8}", fix_lon);
                                self.result = format!(
                                    "DMS Format:\nLatitude: {}\nLongitude: {}\n\nDecimal Format:\nLatitude: {}\nLongitude: {}\n",
                                    dms_lat, dms_lon, dec_lat, dec_lon
                                );
                            },
                            Err(e) => self.error = Some(e),
                        }
                    } else {
                        self.error = Some("Input parsing error".to_string());
                    }
                }

                if let Some(ref err) = self.error {
                    ui.colored_label(egui::Color32::RED, err);
                }

                ui.separator();
                ui.heading("Calculation Results");
                ui.label(&self.result);
            });
        });
    }
}

fn main() {
    let native_options = eframe::NativeOptions::default();
    let _ = eframe::run_native(
        "Fix Calculator",
        native_options,
        Box::new(|_cc| Ok(Box::new(AviationCoordinateCalculatorApp::default()))),
    );
}
