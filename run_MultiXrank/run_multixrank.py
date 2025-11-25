import sys
import time, resource
import multixrank

start = time.time() 

PHENOTYPE = sys.argv[1]

multixrank_obj = multixrank.Multixrank(config=f"{PHENOTYPE}/config.yml", wdir=f"{PHENOTYPE}")

ranking_df = multixrank_obj.random_walk_rank()

multixrank_obj.write_ranking(ranking_df, path=f"{PHENOTYPE}/output_{PHENOTYPE}")

end = time.time()
usage = resource.getrusage(resource.RUSAGE_SELF)
print(f"wall_time={end-start:.6f}s")
print(f"user_time={usage.ru_utime:.6f}s")
print(f"sys_time={usage.ru_stime:.6f}s")
print(f"max_rss_kb={usage.ru_maxrss}")