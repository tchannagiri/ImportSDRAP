def overlaps(a, b):
    if (a['start'] > b['end'] + 1): return False
    if (b['start'] > a['end'] + 1): return False
    return True
  
def get_union(interval_list):
  interval_list = list(sorted(
    interval_list,
    key = lambda a: a['start']
  ))
  interval_list_union = []
  for x in interval_list:
    if (
      (len(interval_list_union) > 0) and
      overlaps(interval_list_union[-1], x)
    ):
      interval_list_union[-1]['end'] = max(interval_list_union[-1]['end'], x['end'])
    else:
      interval_list_union.append(x)
  return interval_list_union