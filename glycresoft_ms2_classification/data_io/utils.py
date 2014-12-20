from multiprocessing import Process, Queue
from Queue import Empty as QueueEmptyException


def _wrap_enqueue_payload(payload_fn, payload_param, out_queue):
    results = payload_fn(*payload_param)
    out_queue.put(results)


class ParallelParser(object):
    def __init__(self, payload_fn, payload_param):
        self.payload_fn = payload_fn
        self.payload_param = payload_param
        self.results_queue = Queue()
        self.worker = Process(target=_wrap_enqueue_payload, args=(payload_fn, payload_param, self.results_queue))
        self.worker.start()
        self.done = False

    def await(self, timeout=None):
        if self.done:
            raise Exception("Process Already Complete")
        while(True):
            try:
                parser = self.results_queue.get()
                break
            except QueueEmptyException:
                pass
        if self.worker.is_alive():
            self.worker.terminate()
        self.done = True
        return parser
