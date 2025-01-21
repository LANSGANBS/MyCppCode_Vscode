fun solve() {
  val n = readLine() !!.toInt() var ans = 0 var allOnes = false var allZeros =
      false var minDiff = 1000000000
      repeat(n) {
    val line = readLine() !!val ones = line.count{it == '1'} val zeros =
        line.length - ones
                          minDiff = minOf(minDiff, abs(ones - zeros)) ans +=
        minOf(ones, zeros)
            when {
      ones > zeros->allOnes = true zeros > ones->allZeros = true else->{
        allOnes = true allZeros = true
      }
    }
  }
  if (ans == 0 || (allZeros && allOnes)) {
    println(ans)
  } else {
    println(ans + minDiff)
  }
}

fun main() {
  val t = 1 repeat(t) { solve() }
}

fun abs(x : Int) : Int = if (x >= 0) x else - x