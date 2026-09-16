from date_utilities import decimal_to_date, date_to_decimal
import pandas as pd


def year_decimal_to_date_tick_labels(series, tick_freq='automatic'):
    """
    Uses year_decimal values to propose date tick labels.

    Parameters
    ----------
    series: pd.Series
        Series of year_decimal values.
    tick_freq: str, default='automatic'
        Suggested tick frequency. Options are 'automatic', 'yearly', 'quarterly',
        'monthly', 'half month' or 'weekly'.

    Returns
    -------
    tick_year_decimals: pd.Series
        Location of ticks. A Series of tick year_decimal values.
    tick_labels:
        The corresponding tick labels.
    """
    max_date = decimal_to_date(series.max())
    min_date = decimal_to_date(series.min())
    if tick_freq == 'automatic':
        time_length = series.max() - series.min()
        if time_length > 5:
            tick_freq = 'yearly'
        elif time_length > 2:
            tick_freq = 'quarterly'
        elif time_length > 1:
            tick_freq = 'monthly'
        elif time_length > 0.5:
            tick_freq = 'half month'
        else:
            tick_freq = 'weekly'

    if tick_freq == 'yearly':
        freq = 'YS'
        strftime = '%Y'
    elif tick_freq == 'quarterly':
        freq = 'QS'
        strftime = '%Y-%b'
    elif tick_freq == 'monthly':
        freq = 'MS'
        strftime = '%Y-%b'
    elif tick_freq == 'half month':
        freq = 'SME'
        strftime = '%Y-%b-%d'
    elif tick_freq == 'weekly':
        freq = 'W'
        strftime = '%Y-%b-%d'
    else:
        raise ValueError \
            ('x_tick_freq must be "automatic", "yearly", "quarterly", "monthly", ''half month" or "weekly".' +
                         '\n x_tick_freq=' + str(tick_freq))

    tick_dates = pd.date_range(start=min_date.replace(day=1), end=max_date, freq=freq)
    tick_labels = tick_dates.to_series().dt.strftime(strftime).to_list()
    tick_year_decimals = tick_dates.map(date_to_decimal)
    return tick_year_decimals, tick_labels