# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.
# 

from gnuradio import gr, gr_unittest
from gnuradio import blocks
import numpy
from time import sleep
import synchronizer_v2
import pmt
import scipy.io as sio
import os

class qa_synchronizer_v2 (gr_unittest.TestCase):

    def setUp (self):
        sps = 2
        Fs = 100000.0
        self.syncBlock = synchronizer_v2.synchronizer_v2(sps, Fs)
        
        print os.getcwd()
        
        x = sio.loadmat('../../matlab/gr_impl_test.mat')
        dt = numpy.dtype(numpy.complex64)
        self.burst1 = x['burst1'].transpose().astype(dt)
        self.burst2 = x['burst2'].transpose().astype(dt)
        self.burst3 = x['burst3'].transpose().astype(dt)
        self.burst4 = x['burst4'].transpose().astype(dt)
        self.burst5 = x['burst5'].transpose().astype(dt)
        
    def tearDown (self):
        None
    
    def test_001_t (self):
        print 'Running Synchronizer Test 1'
        
        pmtMsg = pmt.cons(pmt.PMT_NIL, pmt.to_pmt(self.burst1))
        self.syncBlock.enableDebugMode('/tmp/gr_test1.mat')
        self.syncBlock.handler(pmtMsg)
        
if __name__ == '__main__':
    gr_unittest.run(qa_synchronizer_v2, "qa_synchronizer_v2.xml")
    
