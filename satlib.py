import numpy as np

def calculate_hist_stats(data_list, sample_per_exp, bins=None):
    """
    ### Calculate the mean, standard deviation, and bin centers for histogram bins.
    
    **Parameters:**
    - `data_list` (list or np.ndarray): List containing all experimental data.
    - `sample_per_exp` (int): Number of data points per experiment.
    - `bins` (int, str, or array-like, optional): Bin edges or number of bins (default is 'auto').

    **Returns:**
    - `mean_counts` (np.ndarray): Mean frequency counts for each bin.
    - `std_counts` (np.ndarray): Standard deviation of frequency counts for each bin.
    - `bin_centers` (np.ndarray): Center coordinates of each bin.
    """
    # - `bins` (np.ndarray): Bin edges used for the histogram.
    # 转换为numpy数组
    data = np.asarray(data_list)
    total_samples = data.size
    
    # 计算实验次数并截断数据
    n_exp = total_samples // sample_per_exp
    if total_samples % sample_per_exp != 0:
        print(f"Warning: 截断 {total_samples - n_exp*sample_per_exp} 个无法分组的数据")
        data = data[:n_exp*sample_per_exp]
    
    # 分割为多组实验数据
    data_groups = data.reshape(n_exp, sample_per_exp)
    
    # 确定分箱边界（基于全部数据）
    if bins is None:
        _, bins = np.histogram(data, bins='auto')  # 自动分箱
    else:
        bins = np.asarray(bins)  # 支持自定义分箱
    
    # 计算每组实验的分箱频数
    all_counts = []
    for group in data_groups:
        counts, _ = np.histogram(group, bins=bins)
        all_counts.append(counts)
    all_counts = np.array(all_counts)
    
    # 统计均值和标准差
    mean_counts = np.mean(all_counts, axis=0)
    std_counts = np.std(all_counts, axis=0)
    
    # 计算分箱中心坐标
    bin_centers = (bins[:-1] + bins[1:]) / 2
    # bins = np.asarray(bins)
    return mean_counts, std_counts, bin_centers