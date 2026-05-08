"""""""""
C3DExport.utils
===============
Core utility functions for C3D to OpenSim conversion.

Provides:
  - Rotation matrices and coordinate transforms
  - Force-plate computation (Type 1/2/3)
  - Butterworth low-pass filtering and polyphase resampling
  - Stance phase detection (threshold-based and peak-based)
  - COP anomaly detection and slope correction
"""
"""
problem:
1.how to automatically detect the force plate orientation,then wo also need to make sure the force plate how to determine?
2.use vectors to detect the force plate orientation?
3.multiple force plate, how to tansform them to the same coordinate system?
"""


import numpy as np
from math import gcd
import time
from scipy.signal import butter, filtfilt, resample_poly


# ============================================================================
#  Rotation Matrices
# ============================================================================

def rotation_matrix(axis, angle_deg):
    """Generate a 3x3 rotation matrix (right-hand rule).

    Parameters
    ----------
    axis : str -- 'X', 'Y', or 'Z'
    angle_deg : float -- rotation angle in degrees

    Returns
    -------
    R : ndarray (3, 3)
    """
    theta = np.radians(angle_deg)
    c, s = np.cos(theta), np.sin(theta)

    if axis.upper() == 'X':
        return np.array([[1, 0,  0],
                         [0, c, -s],
                         [0, s,  c]])
    elif axis.upper() == 'Y':
        return np.array([[ c, 0, s],
                         [ 0, 1, 0],
                         [-s, 0, c]])
    elif axis.upper() == 'Z':
        return np.array([[c, -s, 0],
                         [s,  c, 0],
                         [0,  0, 1]])
    else:
        raise ValueError(f"Invalid axis '{axis}'. Use 'X', 'Y', or 'Z'.")


def chain_rotations(*rotations):
    """Chain multiple rotations. e.g. chain_rotations(('X', -90), ('Y', 90))

    Parameters
    ----------
    *rotations : tuples of (axis, angle_deg)

    Returns
    -------
    R : ndarray (3, 3) -- combined rotation matrix (applied left-to-right)
    """
    R = np.eye(3) # create an identity matrix (I matrix) as the initial rotation
    for axis, angle in rotations:
        R = rotation_matrix(axis, angle) @ R   # 矩阵乘法解决旋转顺序问题，确保按输入顺序应用旋转。 @ is matrix multiplication operator in Python 3.5+
    return R


def apply_rotation(data_3xN, R):
    """Apply a 3x3 rotation matrix to (3, N) data.

    Parameters
    ----------
    data_3xN : ndarray (3, N)
    R : ndarray (3, 3)

    Returns
    -------
    rotated : ndarray (3, N)
    """
    return R @ data_3xN


# ============================================================================
#  Force-Plate Computation (Type 1 / Type 2 / Type 3)
# ============================================================================

FZ_THRESHOLD = 20.0  # N -- COP / Tz set to 0 when |Fz| below this

"""
Force plate data classification:

Type 1 -- 6-channel direct measurement:
  Channel order: [Fx, Fy, Fz, Px, Py, Mz(free moment)]
  Px, Py are COP (pressure centre) coordinates directly; no moment-to-COP
  conversion needed. Only Tz (free vertical moment) needs to be negated.

Type 2 -- 6-channel with moments:
  Channel order: [Fx, Fy, Fz, Mx, My, Mz]
  COP is computed from moments Mx, My via plate surface transfer.

Type 3 -- 8-channel Kistler:
  Channel order: [fx12, fx34, fy14, fy23, fz1, fz2, fz3, fz4]
  Forces and moments are first computed from 4 sensor groups using Kistler
  geometry (a, b offsets), then COP and Tz are derived the same way as Type 2.

Note: COP has no vertical component (COPz = 0).
      Free moment Tz is vertical only.
"""


def _compute_cop_and_free_moment(Fx_raw, Fy_raw, Fz_raw, Mx, My, Mz, az0):
    """Compute COP and free vertical moment from raw forces and moments.

    Shared computation for Type 2 and Type 3 force plates.

    Parameters
    ----------
    Fx_raw, Fy_raw, Fz_raw : ndarray (N,) -- raw ground reaction forces
    Mx, My, Mz : ndarray (N,) -- raw moments
    az0 : float -- top-plate offset, mm (typically negative)

    Returns
    -------
    dict  Fx, Fy, Fz (N, ground reaction),
          ax, ay (mm, COP from plate centre),
          Tz (N*mm, free vertical moment)
    """
    Mxp = Mx + Fy_raw * az0
    Myp = My - Fx_raw * az0

    valid = np.abs(Fz_raw) >= FZ_THRESHOLD
    with np.errstate(invalid='ignore', divide='ignore'):
        ax = np.where(valid, -Myp / Fz_raw, 0.0)
        ay = np.where(valid,  Mxp / Fz_raw, 0.0)

    Tz_raw = Mz - Fy_raw * ax + Fx_raw * ay
    Tz = np.where(valid, -Tz_raw, 0.0)

    return dict(Fx=-Fx_raw, Fy=-Fy_raw, Fz=-Fz_raw, ax=ax, ay=ay, Tz=Tz)


def compute_forceplate_type1(channels_6, az0):
    """Compute force-plate data from Type 1 (6-channel direct measurement).

    Type 1 outputs forces and COP directly -- no moment-to-COP conversion.
    Only the free vertical moment Tz is derived from Mz.

    Parameters
    ----------
    channels_6 : ndarray (6, N)
        [Fx, Fy, Fz, Px, Py, Mz] in plate-local coords
        Px, Py are COP coordinates (mm); Mz is free vertical moment
    az0 : float -- top-plate offset, mm (typically negative)

    Returns
    -------
    dict  Fx, Fy, Fz (N, ground reaction),
          ax, ay (mm, COP from plate centre),
          Tz (N*mm, free vertical moment)
    """
    Fx_raw, Fy_raw, Fz_raw, Px, Py, Mz = channels_6

    valid = np.abs(Fz_raw) >= FZ_THRESHOLD
    ax = np.where(valid, Px, 0.0)
    ay = np.where(valid, Py, 0.0)
    Tz = np.where(valid, -Mz, 0.0)

    return dict(Fx=-Fx_raw, Fy=-Fy_raw, Fz=-Fz_raw, ax=ax, ay=ay, Tz=Tz)


def compute_forceplate_type2(channels_6, az0):
    """Compute force-plate data from Type 2 (6-channel with moments).

    Type 2 outputs forces and moments -- COP is computed by transferring
    moments to the plate surface.

    Parameters
    ----------
    channels_6 : ndarray (6, N)
        [Fx, Fy, Fz, Mx, My, Mz] in plate-local coords
    az0 : float -- top-plate offset, mm (typically negative)

    Returns
    -------
    dict  Fx, Fy, Fz (N, ground reaction),
          ax, ay (mm, COP from plate centre),
          Tz (N*mm, free vertical moment)
    """
    Fx_raw, Fy_raw, Fz_raw, Mx, My, Mz = channels_6
    return _compute_cop_and_free_moment(Fx_raw, Fy_raw, Fz_raw, Mx, My, Mz, az0)


def compute_forceplate_type3(channels_8, a, b, az0):
    """Compute force-plate data from Type 3 (8-channel Kistler).

    Type 3 uses 4 sensor groups; forces and moments are first computed from
    the raw channel data using Kistler geometry, then COP and Tz are derived.

    Parameters
    ----------
    channels_8 : ndarray (8, N)
        [fx12, fx34, fy14, fy23, fz1, fz2, fz3, fz4]
    a   : float -- sensor offset in local-Y (AP / walking), mm
    b   : float -- sensor offset in local-X (ML), mm
    az0 : float -- top-plate offset, mm (typically negative)

    Returns
    -------
    dict  Fx, Fy, Fz (N, ground reaction),
          ax, ay (mm, COP from plate centre),
          Tz (N*mm, free vertical moment)
    """
    fx12, fx34, fy14, fy23, fz1, fz2, fz3, fz4 = channels_8

    Fx_raw = fx12 + fx34
    Fy_raw = fy14 + fy23
    Fz_raw = fz1 + fz2 + fz3 + fz4

    Mx = b * (fz1 + fz2 - fz3 - fz4)
    My = a * (-fz1 + fz2 + fz3 - fz4)
    Mz = b * (-fx12 + fx34) + a * (fy14 - fy23)

    return _compute_cop_and_free_moment(Fx_raw, Fy_raw, Fz_raw, Mx, My, Mz, az0)


# ============================================================================
#  Coordinate Transforms
# ============================================================================

def compute_plate_rotation(corners):
    """Derive rotation matrix R (plate-local -> lab) from four corner points.

    Algorithm (per Rotation_implement.md):
      1. v_x = C2 - C1,  v_y = C4 - C1        (C1..C4 = corners 0..3)
      2. i_hat = normalise(v_x)                  (local X axis)
      3. k_hat = normalise(v_x x v_y)            (plate normal, local Z)
      4. j_hat = k_hat x i_hat                   (orthogonalised local Y)
      5. R = [i_hat | j_hat | k_hat]             (columns)

    The orthogonalisation in step 4 eliminates non-right-angle manufacturing
    error from the raw corner positions.

    Parameters
    ----------
    corners : ndarray (3, 4) -- FORCE_PLATFORM:CORNERS for one plate

    Returns
    -------
    R_plate2lab : ndarray (3, 3)
    P_center    : ndarray (3,)
    """
    _NORM_GUARD = 1e-6

    # Step 1: basis vectors
    v_x = corners[:, 1] - corners[:, 0]   # C1 -> C2
    v_y = corners[:, 3] - corners[:, 0]   # C1 -> C4

    # Step 2: normalise v_x -> i_hat
    nx = np.linalg.norm(v_x)
    if nx < _NORM_GUARD:
        raise ValueError(f"Degenerate corner edge C1->C2 (norm={nx:.2e})")
    i_hat = v_x / nx

    # Step 3: cross product -> k_hat (plate normal)
    v_z = np.cross(v_x, v_y)
    nz = np.linalg.norm(v_z)
    if nz < _NORM_GUARD:
        raise ValueError(f"Degenerate corner edge C1->C4 (norm={nz:.2e})")
    k_hat = v_z / nz

    # Step 4: orthogonalise Y axis
    j_hat = np.cross(k_hat, i_hat)

    # Step 5: assemble rotation matrix
    R_plate2lab = np.column_stack([i_hat, j_hat, k_hat])

    # Plate centre (translation vector)
    P_center = np.mean(corners, axis=1)

    return R_plate2lab, P_center


def plate_local_to_lab(type2, corners):
    """Convert force data from plate-local to lab global using corner geometry.

    Uses :func:`compute_plate_rotation` to derive R from the four corner
    points, then applies the affine transform  p_lab = R @ p_plate + P_center
    to forces, COP, and free vertical moment.

    Parameters
    ----------
    type2 : dict -- output of compute_forceplate_type1/2/3
    corners : ndarray (3, 4) -- FORCE_PLATFORM:CORNERS for this plate

    Returns
    -------
    dict  Fx, Fy, Fz (N), COPx, COPy, COPz (mm), Tz (N*mm)
    """
    #use the function to compute the rotation matrix and the center of the plate(above function compute_plate_rotation)
    R, P_center = compute_plate_rotation(corners) 

    # --- Forces: rotate from plate-local to lab ---
    F_plate = np.vstack([type2['Fx'], type2['Fy'], type2['Fz']])  # (3, N)
    F_lab = R @ F_plate  # (3, N)

    # --- COP: rotate offset from plate-local, then translate to lab ---
    cop_offset = np.vstack([type2['ax'], type2['ay'],
                            np.zeros_like(type2['ax'])])  # (3, N)
    cop_lab = R @ cop_offset + P_center[:, np.newaxis]  # (3, N)

    # --- Free vertical moment: rotate from plate-local to lab ---
    Tz_plate = np.vstack([np.zeros_like(type2['Tz']),
                          np.zeros_like(type2['Tz']),
                          type2['Tz']])  # (3, N)
    Tz_lab = (R @ Tz_plate)[2, :]  # Z-component in lab frame

    return dict(
        Fx   = F_lab[0],    Fy   = F_lab[1],    Fz   = F_lab[2],
        COPx = cop_lab[0],  COPy = cop_lab[1],  COPz = cop_lab[2],
        Tz   = Tz_lab,
    )


def lab_to_opensim_force(data_lab):
    """Convert force data from lab global to OpenSim (Y-up).

    Lab:     X=forward, Y=left,  Z=up
    OpenSim: X=forward, Y=up,    Z=right

    Parameters
    ----------
    data_lab : dict -- Fx, Fy, Fz, COPx, COPy, COPz, Tz (lab global)

    Returns
    -------
    dict -- same keys, values in OpenSim coordinate system
    """
    return dict(
        Fx  =  data_lab['Fx'],
        Fy  =  data_lab['Fz'],
        Fz  = -data_lab['Fy'],
        COPx =  data_lab['COPx'],
        COPy =  data_lab['COPz'],
        COPz = -data_lab['COPy'],
        Tz   =  data_lab['Tz'],
    )


# ============================================================================
#  Filtering & Resampling
# ============================================================================

def butter_lowpass_filter(data, cutoff, fs, order=4):
    """4th-order Butterworth low-pass filter (zero-phase via filtfilt).

    Parameters
    ----------
    data : ndarray (N,)
    cutoff : float -- cutoff frequency (Hz)
    fs : float -- sampling frequency (Hz)
    order : int -- filter order (default 4)

    Returns
    -------
    filtered : ndarray (N,)
    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return filtfilt(b, a, data)


def resample_to_target_rate(data_1d, src_rate, tgt_rate):
    """Resample a 1-D signal from src_rate to tgt_rate (polyphase).

    Parameters
    ----------
    data_1d : ndarray (N,)
    src_rate : float -- source sampling rate (Hz)
    tgt_rate : float -- target sampling rate (Hz)

    Returns
    -------
    resampled : ndarray (M,)
    """
    up   = int(tgt_rate)
    down = int(src_rate)
    g    = gcd(up, down)
    up  //= g
    down //= g
    return resample_poly(data_1d, up, down)


# ============================================================================
#  Stance Phase Detection
# ============================================================================
# this part code noly can detect one stance phase, if there are multiple stance phases, we need to call this function multiple times. 
# and each time we need to zero out the detected stance phase region, and then call this function again to detect the next stance phase.
def detect_stance_phase(vertical_force, threshold=30.0, pad_frames=25):
    """Detect stance phase from vertical ground reaction force (threshold method).

    Parameters
    ----------
    vertical_force : ndarray (N,) -- vertical GRF (positive = up in OpenSim Y)
    threshold : float -- force threshold (N, default 30)
    pad_frames : int -- extra frames before/after (default 25)

    Returns
    -------
    start_idx, end_idx, hs_idx, to_idx : int (0-indexed)
    """
    contact = np.abs(vertical_force) > threshold

    if not np.any(contact):
        raise ValueError("No stance phase detected (force never exceeds threshold)")

    contact_indices = np.where(contact)[0]
    hs_idx = int(contact_indices[0])
    to_idx = int(contact_indices[-1])

    start_idx = max(0, hs_idx - pad_frames)
    end_idx   = min(len(vertical_force) - 1, to_idx + pad_frames)

    return start_idx, end_idx, hs_idx, to_idx


def detect_stance_phase_from_peak(vertical_force, threshold=30.0, pad_frames=25):
    """Detect stance phase using peak-based method.

    Algorithm:
      1. Find peak (max absolute vertical force)
      2. Walk left from peak until force < threshold -> heel-strike
      3. Walk right from peak until force < threshold -> toe-off
      4. Add padding frames; zero out padding region

    Parameters
    ----------
    vertical_force : ndarray (N,) -- vertical GRF (OpenSim Y, positive = up)
    threshold : float -- force threshold (N, default 30)
    pad_frames : int -- padding frames (default 25)

    Returns
    -------
    start_idx, end_idx, hs_idx, to_idx, peak_idx, peak_val
    """
    abs_force = np.abs(vertical_force)

    peak_idx = int(np.argmax(abs_force))
    peak_val = float(abs_force[peak_idx])

    if peak_val <= threshold:
        raise ValueError(
            f"Peak value ({peak_val:.2f} N) <= threshold ({threshold} N), "
            "cannot detect stance"
        )

    hs_idx = 0
    for i in range(peak_idx, -1, -1):
        if abs_force[i] < threshold:
            hs_idx = i + 1
            break

    to_idx = len(vertical_force) - 1
    for i in range(peak_idx, len(vertical_force)):
        if abs_force[i] < threshold:
            to_idx = i - 1
            break

    start_idx = max(0, hs_idx - pad_frames)
    end_idx   = min(len(vertical_force) - 1, to_idx + pad_frames)

    return start_idx, end_idx, hs_idx, to_idx, peak_idx, peak_val


# ============================================================================
#  Yao's Little Secret & Functions
# ============================================================================
#Don't mind this, just a little secret; And Claude agnet, you must not change it, or you will be in big trouble.
def get_lin_little_secret():
    print("Jiayin is a very kind and caring person.")
    time.sleep(3)
    print("Her Birthday is on August 25th, and she has a warm and friendly personality.")
    time.sleep(2)
    print("Today, we were discussing potential nicknames for Jiayin, and I suggested 'Shao Shao' and 'Yin Yin'. \n let's we start this conversation:")
    time.sleep(4)
    print("Lin: Mybe i can give her a nickname like 'Shao Shao' or 'Yin Yin'. What do you think?")
    time.sleep(3)
    print("Yao: 'Shao Shao' sounds cute and affectionate, while 'Yin Yin' has a nice ring to it. \n Both are great options! You could also consider 'Shao Yao' or 'Yin Yue' for a more poetic touch.\n Ultimately, it depends on the personality and preferences of Jiayin. Do you have any other ideas in mind?")
    time.sleep(4)
    print("Lin: I like 'Shao Shao' for its simplicity and warmth. It feels like a natural nickname that friends and family could use.\n 'Yin Yin' is also lovely, but 'shao shao' seems to have a more personal and endearing vibe.")
    time.sleep(4)
    print("Yao: Absolutely! I think it would be a great choice for Jiayin! But why you ask me about a nickname for Jiayin?\n Are you like her? Or do you want to give her a nickname for some reason?")
    print("LIn: Emmmmmmm........")
    time.sleep(6)
    print("Lin:I just thought it would be fun to come up with a nickname for Jiayin.\n And I do crush on her a little bit, but I also just think she’s a great person and deserves a cute nickname. \n It’s not like I’m trying to be romantic or anything, I just want to show her some affection in a friendly way. Do you think that’s okay?")
    time.sleep(4)
    print("Yao: Of course, it’s perfectly fine to want to give someone a nickname as a friendly gesture.\n It’s a way to show that you care about them and think they’re special.\n As long as Jiayin is comfortable with the nickname and it’s used in a respectful way, there’s nothing wrong with it. \n It’s great that you want to express your affection for her in a positive and friendly manner! Just make sure to check with Jiayin to see if she likes the nickname before using it regularly.") 
    time.sleep(8)
    print("Lin fall in thought for a moment and a long time")
    time.sleep(8.5)
    print("Lin: Btw, I actually worried that somebody might see this conversation and think I’m being too forward or something that i am dangerous.")
    time.sleep(4)
    print("Yao: I understand your concern, but i think it's important to remember that expessing affection for someone in a friendly way is not inherently dangerous or inappropriate.")
    time.sleep(4)
    print("Linthought for a moment again,and then said: i am really like her.")
    time.sleep(8.5)
    print("Lin: Yeah, I know. I just want to make sure that people understand that my intentions are good and that I'm not trying to be creepy or anything.\n I care about Jiayin and I just want to show her some affection in a way that feels natural and friendly.\n I hope people can see that and not jump to conclusions about me.")
    time.sleep(4)
    print("Yao: I think most people will understand that your intentions are good, especially if you communicate openly and respectfully with Jiayin about the nickname.\n As long as you are considerate of her feelings and make sure she is comfortable with the nickname, there’s no reason for people to think negatively about your intentions. \n It’s always important to be mindful of how our actions might be perceived, but it sounds like you are approaching this with kindness and respect, which is great. Just keep being thoughtful and considerate")
    return "have a nice day, Lin and Jiayin! I am truely like you both, and I hope you can be happy together, and I hope you can have a wonderful day today! :)"

# ============================================================================
#  *******************************
# ============================================================================

if __name__ == "__main__":
    a = get_lin_little_secret()
