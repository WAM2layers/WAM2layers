from pandas import Timestamp, Timedelta

def timeloop(
    start: Timestamp,
    end: Timestamp,
    step: Timedelta,
):
    """Iterator returning tuples of (t, t+1/2, t+1).

    Example:

        >>> from pandas import Timestamp, Timedelta
        >>>
        >>> start = Timestamp(2023, 1, 1)
        >>> end = Timestamp('20230102T1800')
        >>> freq = Timedelta('6h')
        >>>
        >>> time = timeloop(start, end, freq)
        >>> for t, th, tn in time:
        ...     print(f"{t}, {th}, {tn}")
        2023-01-01 00:00:00, 2023-01-01 03:00:00, 2023-01-01 06:00:00
        2023-01-01 06:00:00, 2023-01-01 09:00:00, 2023-01-01 12:00:00
        2023-01-01 12:00:00, 2023-01-01 15:00:00, 2023-01-01 18:00:00
        2023-01-01 18:00:00, 2023-01-01 21:00:00, 2023-01-02 00:00:00
        2023-01-02 00:00:00, 2023-01-02 03:00:00, 2023-01-02 06:00:00
        2023-01-02 06:00:00, 2023-01-02 09:00:00, 2023-01-02 12:00:00
        2023-01-02 12:00:00, 2023-01-02 15:00:00, 2023-01-02 18:00:00

    """
    if start > end:
        raise ValueError("end should be later than start")

    t = start
    while t < end:
        yield t, t + step / 2, t + step
        t += step


def reverse_timeloop(
    start: Timestamp,
    end: Timestamp,
    step: Timedelta,
):
    """Iterator returning tuples of (t, t-1/2, t-1).

    Example:

        >>> from pandas import Timestamp, Timedelta
        >>>
        >>> start = Timestamp(2023, 1, 1)
        >>> end = Timestamp('20230102T1800')
        >>> freq = Timedelta('6h')
        >>>
        >>> time = reverse_timeloop(start, end, freq)
        >>> for t, th, tn in time:
        ...     print(f"{t}, {th}, {tn}")
        2023-01-02 18:00:00, 2023-01-02 15:00:00, 2023-01-02 12:00:00
        2023-01-02 12:00:00, 2023-01-02 09:00:00, 2023-01-02 06:00:00
        2023-01-02 06:00:00, 2023-01-02 03:00:00, 2023-01-02 00:00:00
        2023-01-02 00:00:00, 2023-01-01 21:00:00, 2023-01-01 18:00:00
        2023-01-01 18:00:00, 2023-01-01 15:00:00, 2023-01-01 12:00:00
        2023-01-01 12:00:00, 2023-01-01 09:00:00, 2023-01-01 06:00:00
        2023-01-01 06:00:00, 2023-01-01 03:00:00, 2023-01-01 00:00:00

    """
    if start > end:
        raise ValueError("end should be later than start")

    t = end
    while t > start:
        yield t, t - step / 2, t - step
        t -= step


if __name__=="__main__":

    start = Timestamp(2023, 1, 1)
    end = Timestamp('20230102T1800')
    freq = Timedelta('6h')

    print('Timeloop')
    time = timeloop(start, end, freq)
    for t, th, tn in time:
        print(f"{t}, {th}, {tn}")

    print('\nReverse timeloop')
    time = reverse_timeloop(start, end, freq)
    for t, th, tn in time:
        print(f"{t}, {th}, {tn}")
