# Снимок клиентской области окна molview плюс увеличенный фрагмент.
# Вспомогательный сценарий для разработки, в проект не входит.
param(
  [string]$Shot = "snap.png",
  [int]$ZoomX = -1, [int]$ZoomY = 0, [int]$ZoomW = 340, [int]$ZoomH = 260, [int]$Factor = 3
)

Add-Type -AssemblyName System.Drawing
$src = @'
using System;
using System.Runtime.InteropServices;
public class Snp {
  [DllImport("user32.dll")] public static extern bool GetClientRect(IntPtr h, out RECT r);
  [DllImport("user32.dll")] public static extern bool ClientToScreen(IntPtr h, ref POINT p);
  [DllImport("user32.dll")] public static extern bool SetForegroundWindow(IntPtr h);
  [DllImport("user32.dll")] public static extern bool SetWindowPos(IntPtr h, IntPtr a, int x, int y, int cx, int cy, uint f);
  [StructLayout(LayoutKind.Sequential)] public struct RECT { public int Left, Top, Right, Bottom; }
  [StructLayout(LayoutKind.Sequential)] public struct POINT { public int X, Y; }
}
'@
if (-not ("Snp" -as [type])) { Add-Type -TypeDefinition $src }

$p = Get-Process molview
$null = [Snp]::SetWindowPos($p.MainWindowHandle, [IntPtr](-1), 0, 0, 0, 0, 0x0003)
$null = (New-Object -ComObject WScript.Shell).AppActivate($p.Id)
$null = [Snp]::SetForegroundWindow($p.MainWindowHandle)
Start-Sleep -Milliseconds 900

$r = New-Object Snp+RECT; $null = [Snp]::GetClientRect($p.MainWindowHandle, [ref]$r)
$o = New-Object Snp+POINT; $null = [Snp]::ClientToScreen($p.MainWindowHandle, [ref]$o)
$bmp = New-Object System.Drawing.Bitmap $r.Right, $r.Bottom
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.CopyFromScreen($o.X, $o.Y, 0, 0, $bmp.Size)
$bmp.Save("C:\git\fun\molviewcpp\build\$Shot", [System.Drawing.Imaging.ImageFormat]::Png)

if ($ZoomX -ge 0) {
  $crop = New-Object System.Drawing.Rectangle $ZoomX, $ZoomY, $ZoomW, $ZoomH
  $part = $bmp.Clone($crop, $bmp.PixelFormat)
  $zoom = New-Object System.Drawing.Bitmap ($ZoomW * $Factor), ($ZoomH * $Factor)
  $gz = [System.Drawing.Graphics]::FromImage($zoom)
  $gz.InterpolationMode = [System.Drawing.Drawing2D.InterpolationMode]::NearestNeighbor
  $gz.PixelOffsetMode = [System.Drawing.Drawing2D.PixelOffsetMode]::Half
  $gz.DrawImage($part, 0, 0, $ZoomW * $Factor, $ZoomH * $Factor)
  $zoom.Save("C:\git\fun\molviewcpp\build\zoom_$Shot", [System.Drawing.Imaging.ImageFormat]::Png)
}
"снимок $Shot готов"
